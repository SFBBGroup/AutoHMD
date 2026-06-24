# AutoHMD — Automated Heated Molecular Dynamics Pipeline for Antibody-Antigen Complexes

AutoHMD is a modular pipeline for the automated preparation and execution of Heated Molecular Dynamics (HMD) simulations using the AMBER suite. Optimized for protein–antibody systems, it streamlines the entire workflow—from raw PDB input to sequential multi-temperature production (310 K → 330 K → 360 K → 390 K) and advanced trajectory analysis.

The pipeline automates labor-intensive tasks such as system protonation, terminal capping, and Hydrogen Mass Repartitioning (HMR), while facilitating comprehensive post-processing, including iRMSD calculations, automated plotting, and optional interaction profiling via PLIP integration.

The curated antibody-antigen structure and HMD trajectories dataset, containing over 300 simulations, is available on Zenodo.

---
## Citation

If you use AutoHMD in your work, please cite:

> *(citation placeholder — update before publication)*


---


## Requirements

AutoHMD was implemented and tested using **Python 3.12**. The following external softwares are required:

| Software | Tested Version | Role |
|---|---|---|
| AMBER / AmberTools | [Amber24](https://ambermd.org/AmberMD.php) | MD engine, `tleap`, `parmed`, `cpptraj`, `pmemd.cuda` |
| PDB2PQR | [3.6.1](https://github.com/Electrostatics/pdb2pqr/releases/tag/v3.6.1) | Protonation state assignment via PROPKA |
| PDBTools | [2.5.1](https://github.com/haddocking/pdb-tools) | PDB manipulation (`pdb_reres`, `pdb_fixinsert`) |
| PLIP | [2.3.1](https://github.com/pharmai/plip/releases/tag/v2.3.1) | Protein–ligand interaction analysis *(optional, additional_analysis only)* |

By default, AutoHMD loads these dependencies as **environment modules** (`module load`). Users may alternatively modify the execution scripts to use locally installed software, provided all required executables are accessible via `PATH`.

Python packages required for the auxiliary and analysis scripts:

```
numpy (tested in v2.3.0)
matplotlib (tested in v3.10.7)
biopython
pyyaml (tested in v6.0.3)
pandas
seaborn
tqdm (tested in v4.67.1)
```

---

## Quick Start

### 1. Clone the repository

```bash
git clone https://github.com/SFBBGroup/AutoHMD.git
cd AutoHMD
```

### 2. Place your PDB file in `prepare/`

```bash
cp /path/to/your_protein.pdb prepare/
```

### 3. Edit `prepare/config_hMD.conf`

The only required change is the protein filename. 

```bash
INPUT_PDB="your_protein"   # filename without .pdb extension
FORCEFIELD="ff19SB"
WATER_MODEL="opc"
PH="7.4"
SOLVATION_DISTANCE="12.0"
ION_CONCENTRATION="0.15"
```

### 4. Run system preparation

```bash
cd prepare/
bash prepare_autohMD.sh
```

On success, the files `HMR_solvated_system.prmtop`, `solvated_system.rst7`, and `solvated_system.pdb` are automatically copied to `submit/`.

### 5. Submit or run the MD simulation

```bash
cd ../submit/

# Via HTCondor (3 replicas in parallel):
condor_submit hMD_jobSubmit

# Or locally (single run):
bash executable.sh
```

### 6. (Optional) Run PLIP interaction analysis

```bash
cd ../additional_analysis/
# Edit config.yaml with your trajectory paths, then:
bash submit_plip_analysis.sh
```


---


