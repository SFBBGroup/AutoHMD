# submit/ — Simulation Submission and Execution

This folder contains the files needed to run the full Heated MD simulation protocol, from minimization through multi-temperature production MD, followed by automatic trajectory compaction, iRMSD calculation, and RMSD plotting.

By default, `executable.sh` loads modules:

```bash
module load amber/A24
module load python/3.12.3/gcc13.2/miniconda_24.4.0
```

Modify these lines if your environment uses different module names or locally installed software.

## Usage

### HTCondor submission

The provided `hMD_jobSubmit` submits **3 independent replicas** in parallel, each on a separate GPU node. Outputs from each replica are collected in `rep0/`, `rep1/`, and `rep2/`.

```bash
condor_submit hMD_jobSubmit
```

Key settings in `hMD_jobSubmit`:

| Setting | Value | Description |
|---|---|---|
| `queue` | `3` | Number of replicas to submit |
| `request_gpus` | `1` | GPUs per job |
| `request_memory` | `4G` | RAM per job |
| `transfer_input_files` | see file | Files sent to each worker node |
| `transfer_output_remaps` | `rep=rep$(Process)` | Renames output folder per replica |

To restrict jobs to a specific machine, uncomment and edit:

```
requirements = (Machine == "your_cluster")
```

### Local execution

To run a single replica locally (no HTCondor):

```bash
bash executable.sh
```

Output will be written to `rep/` in the current directory.

---

## Simulation Protocol

`executable.sh` creates all AMBER input (`.in`) files at runtime and runs the following stages sequentially using `pmemd.cuda`:

| Step | File | Type | Length | Temperature | Restraints |
|---|---|---|---|---|---|
| 1a | `1_water_minimization.in` | Minimization | 10 000 cycles | — | Solute fixed (1 kcal/mol/Å²) |
| 1b | `1_part2_system_minimization.in` | Minimization | 20 000 cycles | — | None |
| 2 | `2_heat.in` | NVT heating | 500 ps | 0 → 310 K | Heavy atoms 100 kcal/mol/Å² |
| 3 | `3_Water_density.in` | NPT | 1 ns | 310 K | Heavy atoms 100 kcal/mol/Å² |
| 4 | `4_Water_densitY_soft.in` | NPT | 1 ns | 310 K | Heavy atoms 10 kcal/mol/Å² |
| 5 | `5_2nd_Minimization.in` | Minimization | 10 000 cycles | — | Backbone 10 kcal/mol/Å² |
| 6 | `6_First_equil_10kcalmolA2.in` | NPT | 1 ns | 310 K | Backbone 10 kcal/mol/Å² |
| 7 | `7_Second_equilibration_1kcalmolA2.in` | NPT | 1 ns | 310 K | Backbone 1 kcal/mol/Å² |
| 8 | `8_Third_equilibration_01kcalmolA2.in` | NPT | 1 ns | 310 K | Backbone 0.1 kcal/mol/Å² |
| 9 | `9_Unrestr_equilibration.in` | NPT | 1 ns | 310 K | None |
| 10 | `mdHMR310.in` | NPT production | 30 ns | 310 K | None |
| 11 | `mdHMR330.in` | NPT production | 12.5 ns | 330 K | None |
| 12 | `mdHMR360.in` | NPT production | 12.5 ns | 360 K | None |
| 13 | `mdHMR390.in` | NPT production | 15 ns | 390 K | None |

All production stages use a **4 fs timestep** enabled by the HMR topology. SHAKE constraints are applied to bonds involving hydrogen (`ntc=2`, `ntf=2`). Langevin thermostat (`ntt=3`, `gamma_ln=1.0`) and Berendsen barostat (`ntp=1`, `barostat=1`) are used throughout NPT stages.

---

## Post-Processing

After production MD completes, `executable.sh` automatically runs the following:

### compactTraj_CalcRmsd.py

Strips water and ions from the production trajectories and calculates the interface RMSD (iRMSD) between antibody and antigen using `cpptraj`.

**Fragment detection logic:**

- Fragments starting with `ACE` → **antigen**
- First two fragments not starting with `ACE` → **antibody** (heavy chain H + light chain L)

The iRMSD is calculated on Cα atoms at the interface, defined by a distance cutoff of 8.0 Å between antibody and antigen residues.

Called by `executable.sh` as:

```bash
python3 compactTraj_CalcRmsd.py solvated_system.pdb \
    --prmtop-in HMR_solvated_system.prmtop \
    --traj-in rep/trajectory/mdProd*.nc \
    --prmtop-out rep/analysis/rep0_noWat_hMD.prmtop \
    --traj-out rep/analysis/rep0_noWat_hMD.nc \
    --reference nowat_system.rst7 \
    --save rep/analysis/cpptraj_full.in \
    --output rep/analysis/iRMSD_r0.dat \
    --run
```

Raw trajectory files (`.nc`) are deleted after compaction to save disk space.

**Standalone options:**

```bash
# Show detected fragments and residue masks only:
python3 compactTraj_CalcRmsd.py system.pdb --info

# Generate cpptraj script without executing:
python3 compactTraj_CalcRmsd.py system.pdb [options] --save script.in

# Compact only (skip RMSD):
python3 compactTraj_CalcRmsd.py system.pdb [options] --compact-only

# RMSD only (trajectory already compacted):
python3 compactTraj_CalcRmsd.py system.pdb [options] --rmsd-only
```

---

### RMSD_plot.py

Generates an iRMSD plot for a **single replica** from `rep/analysis/iRMSD_r0.dat`.

```bash
python3 RMSD_plot.py
```

The plot spans 70 ns total with vertical dashed lines marking the temperature transitions (310 K → 330 K → 360 K → 390 K) and a horizontal dashed line at 5 Å. Output is saved to `rep/analysis/iRMSD.png` at 400 dpi.

---

### RMSD_plot_3replicates.py

Generates an iRMSD plot overlaying **three replicas** from `rep0/`, `rep1/`, and `rep2/`.

```bash
python3 RMSD_plot_3replicates.py
```

Expected input files:

```
rep0/analysis/iRMSD_r0.dat
rep1/analysis/iRMSD_r0.dat
rep2/analysis/iRMSD_r0.dat
```

Output is saved to `iRMSD_replicas.png` at 400 dpi. Run this script after all three HTCondor replicas have completed.

---

