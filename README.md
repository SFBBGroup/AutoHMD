# AutoHMD
# AutoHMD — Automated Heated Molecular Dynamics Pipeline

AutoHMD is a modular pipeline for automated preparation and execution of Heated Molecular Dynamics (hMD) simulations using AMBER. It is designed for protein–antibody systems and handles all steps from raw PDB input through production MD, trajectory compaction, RMSD calculation, and optional protein–protein interaction analysis via PLIP.

---

## Table of Contents

- [Overview](#overview)
- [Pipeline Structure](#pipeline-structure)
- [Requirements](#requirements)
- [Quick Start](#quick-start)
- [Documentation](#documentation)
- [Citation](#citation)
- [License](#license)

---

## Overview

AutoHMD runs heated MD simulations by sequentially raising the simulation temperature (310 K → 330 K → 360 K → 390 K), enabling enhanced conformational sampling relative to conventional MD. The pipeline automates:

- System preparation (protonation, gap handling, terminal capping, solvation, ion placement, hydrogen mass repartitioning)
- Simulation execution (minimisation, heating, equilibration, multi-temperature production)
- Post-processing (trajectory compaction, iRMSD calculation, RMSD plotting)
- Optional interaction analysis (PLIP-based protein–protein interaction profiling along the trajectory)

---

## Pipeline Structure

```
AutoHMD/
├── README.md                         ← You are here
│
├── prepare/                          # System preparation
│   ├── README.md
│   ├── config_hMD.conf
│   ├── prepare_autohMD.sh
│   └── scripts/
│       ├── README.md
│       ├── fix_gaps_with_ter.py
│       ├── pdb4amber_ter_fix.py
│       ├── capping.py
│       ├── detect_and_fix_gaps.py
│       └── ion_concentration.py
│
├── submit/                           # Simulation submission and execution
│   ├── README.md
│   ├── hMD_jobSubmit
│   ├── executable.sh
│   ├── compactTraj_CalcRmsd.py
│   ├── RMSD_plot.py
│   └── RMSD_plot_3replicates.py
│
└── additional_analysis/              # Optional post-simulation analysis
    ├── README.md
    ├── config.yaml
    ├── submit_plip_analysis.sh
    └── plip_trajectory_analysis.py
```

---

## Requirements

AutoHMD was implemented and tested using **Python 3.12**. The following external software is required:

| Software | Tested Version | Role |
|---|---|---|
| AMBER / AmberTools | Amber24 | MD engine, `tleap`, `parmed`, `cpptraj`, `pmemd.cuda` |
| PDB2PQR | 3.6.1 | Protonation state assignment via PROPKA |
| PDBTools | 2.5.1 | PDB manipulation (`pdb_reres`, `pdb_fixinsert`) |
| PLIP | 3.12 | Protein–ligand interaction analysis *(optional, additional_analysis only)* |

By default, AutoHMD loads these dependencies as **environment modules** (`module load`). Users may alternatively modify the execution scripts to use locally installed software, provided all required executables are accessible via `PATH`.

Python packages required for the auxiliary and analysis scripts:

```
numpy
matplotlib
biopython
pyyaml
pandas
seaborn
tqdm
```

---

## Quick Start

### 1. Clone the repository

```bash
git clone https://github.com/your-org/AutoHMD.git
cd AutoHMD
```

### 2. Place your PDB file in `prepare/`

```bash
cp /path/to/your_protein.pdb prepare/
```

### 3. Edit `prepare/config_hMD.conf`

The only required change is the protein filename. All other parameters have validated defaults.

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

## Documentation

Each pipeline stage has its own detailed README:

| Folder | README | Contents |
|---|---|---|
| `prepare/` | [`prepare/README.md`](prepare/README.md) | Configuration, preparation steps, output files |
| `prepare/scripts/` | [`prepare/scripts/README.md`](prepare/scripts/README.md) | Individual script usage and arguments |
| `submit/` | [`submit/README.md`](submit/README.md) | Simulation stages, HTCondor setup, post-processing |
| `additional_analysis/` | [`additional_analysis/README.md`](additional_analysis/README.md) | PLIP analysis, chain assignment modes, output |

---

## Citation

If you use AutoHMD in your work, please cite:

> *(citation placeholder — update before publication)*

---

## License

*(license placeholder)*
