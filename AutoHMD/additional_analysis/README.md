# additional_analysis/ — PLIP Trajectory Analysis

This folder contains scripts for protein–protein interaction profiling along the compacted MD trajectories using [PLIP](https://github.com/pharmai/plip) (Protein–Ligand Interaction Profiler). The analysis extracts individual frames from the trajectory, runs PLIP on each frame, and produces interaction frequency tables and timeline plots.


---

## Prerequisites

- Python 3.12
- PLIP 2.3.1
- AmberTools (for `cpptraj`, used to extract frames)
- PDBTools 2.5.1 (for `pdb_reres`, `pdb_fixinsert`)
- A conda environment named `plip_analysis` with the following packages:

```
pyyaml
pandas
matplotlib
seaborn
numpy
tqdm
```

If the `plip_analysis` conda environment does not exist, `submit_plip_analysis.sh` will create it automatically.

By default, the submission script loads:

```bash
module load python/3.12.3/gcc13.2/miniconda_24.4.0
module load pdb-tools/2.5.1
module load PLIP/3.12
module load amber/A24
```

Modify these lines if your environment uses different module names.

---

## Usage

### Config file mode

Edit `config.yaml` with your topology, trajectory paths, chain assignments, and output settings, then run:

```bash
# Using default config.yaml:
bash submit_plip_analysis.sh

# Using a custom config file:
bash submit_plip_analysis.sh my_config.yaml
```

### CLI mode

All parameters can be provided directly on the command line. In CLI mode, chain assignment is always performed automatically (auto-detect).

```bash
bash submit_plip_analysis.sh --cli \
    --topology ../submit/rep0/analysis/rep0_noWat_hMD.prmtop \
    --trajectories ../submit/rep0/analysis/rep0_noWat_hMD.nc \
                   ../submit/rep1/analysis/rep0_noWat_hMD.nc \
                   ../submit/rep2/analysis/rep0_noWat_hMD.nc \
    --interval 32 \
    --output plip_results \
    --replica-size 175
```

| Option | Default | Description |
|---|---|---|
| `--topology` | *(required)* | Topology file (`.prmtop` or `.pdb`) |
| `--trajectories` | *(required)* | One or more trajectory files |
| `--output` | `plip_results` | Output directory |
| `--prefix` | `analysis` | Output file prefix |
| `--interval` | `10` | Extract every Nth frame for analysis |
| `--start-frame` | `1` | First frame to analyse |
| `--end-frame` | `-1` | Last frame to analyse (`-1` = last frame) |
| `--target-chain` | `B` | Chain to analyse interactions against (antigen) |
| `--replica-size` | `175` | Frames per replica, used for timeline plot boundary lines |

---

## Configuration — `config.yaml`

```yaml
# Trajectory files
trajectory:
  topology: "../submit/rep0/analysis/rep0_noWat_hMD.prmtop"
  trajectories:
    - "../submit/rep0/analysis/rep0_noWat_hMD.nc"
    - "../submit/rep1/analysis/rep0_noWat_hMD.nc"
    - "../submit/rep2/analysis/rep0_noWat_hMD.nc"
  reference: null

# Frame selection
frames:
  start: 1
  end: -1
  interval: 32

# Chain assignment (see Chain Assignment Modes below)
auto_detect_chains: true
reference_for_chains: null

# PLIP settings
plip:
  target_chain: "B"

# Output settings
output:
  directory: "plip_trajectory_results"
  prefix: "antibody_antigen"

# Chains to keep in processed PDB frames
chains_to_keep:
  - H
  - L
  - B

# Residue renumbering start per chain
reres_starts:
  H: 1
  L: 1
  B: 1

# Chains requiring pdb_fixinsert
fixinsert_chains:
  - H
  - L

# Frames per replica (for timeline plot vertical lines)
replica_size: 175
```

---

## Chain Assignment Modes

The analysis requires residue ranges to be assigned to the antigen (chain B), heavy chain (H), and light chain (L). Three modes are available:

### Option A — Auto-detect *(recommended)*

```yaml
auto_detect_chains: true
```

Requires that ACE/NME terminal capping was applied during system preparation (default in AutoHMD).

**Detection logic:**

- Fragments starting with `ACE` → **Antigen (chain B)**
- First non-ACE fragment → **Heavy chain (H)**
- Second non-ACE fragment → **Light chain (L)**

### Option B — Manual residue ranges

```yaml
auto_detect_chains: false
residue_ranges:
  B: [1, 566]
  H: [567, 686]
  L: [687, 799]
```

Use this option if the system was not prepared with ACE/NME caps, or if auto-detection gives incorrect assignments.

### Option C — Reference PDB

```yaml
auto_detect_chains: false
reference_for_chains: "reference_with_chains.pdb"
```

Chain IDs are read from a reference PDB file. Useful when the compacted trajectory system was built from a PDB with explicit chain labels.

---

## Output

All output is written to the directory specified by `output.directory` (default: `plip_trajectory_results/`).

| File/Folder | Description |
|---|---|
| `frames/` | Extracted and processed PDB files for each analysed frame |
| `plip_output/` | Raw PLIP output (XML and text reports) per frame |
| `*_interactions.csv` | Tabular summary of all interactions across all frames |
| `*_frequency.csv` | Interaction frequency per residue pair across frames |
| `*_timeline.png` | Heatmap of interaction presence over simulation time |
| `*_frequency_plot.png` | Bar chart of the most frequent interactions |

The timeline plot includes vertical lines at replica boundaries (controlled by `replica_size`) to facilitate comparison across replicas.
