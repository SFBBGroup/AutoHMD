#!/bin/bash
# ============================================================
# submit_plip_analysis.sh — Unified PLIP Trajectory Analysis
# ============================================================
#
# USAGE — TWO MODES:
#
#  Mode 1: Config file (auto-detect OR manual, controlled by config)
#    ./submit_plip_analysis.sh [config.yaml]
#    Default config file: config.yaml
#
#  Mode 2: Full command line (always auto-detect)
#    ./submit_plip_analysis.sh --cli \
#        --topology ../submit/rep0/analysis/rep0_noWat_hMD.prmtop \
#        --trajectories ../rep0/rep0.nc ../rep1/rep1.nc ../rep2/rep2.nc \
#        --interval 32 \
#        --output plip_results \
#        --replica-size 175
#
# ============================================================

set -euo pipefail

# ---- Helper: print usage ----
usage() {
    echo ""
    echo "Usage:"
    echo "  $0                          # uses config.yaml"
    echo "  $0 my_config.yaml           # uses a specific config file"
    echo "  $0 --cli [options...]       # full CLI mode (auto-detect chains)"
    echo ""
    echo "CLI mode options:"
    echo "  --topology FILE             Topology file (.prmtop, .pdb, ...)"
    echo "  --trajectories FILE [...]   One or more trajectory files"
    echo "  --output DIR                Output directory (default: plip_results)"
    echo "  --prefix STR                Output file prefix (default: analysis)"
    echo "  --interval N                Extract every Nth frame (default: 10)"
    echo "  --start-frame N             First frame (default: 1)"
    echo "  --end-frame N               Last frame, -1 = last (default: -1)"
    echo "  --target-chain CHAIN        Chain to analyse against (default: B)"
    echo "  --replica-size N            Frames per replica for plots (default: 175)"
    echo ""
}

# ============================================================
# Parse first argument to decide mode
# ============================================================
CLI_MODE=false
CONFIG_FILE="config.yaml"

if [[ $# -eq 0 ]]; then
    # No arguments → config file mode with default config.yaml
    CLI_MODE=false
elif [[ "$1" == "--cli" ]]; then
    CLI_MODE=true
    shift  # remove --cli; the rest are passed straight to Python
elif [[ "$1" == "--help" || "$1" == "-h" ]]; then
    usage
    exit 0
elif [[ "$1" != --* ]]; then
    # Positional argument → treat as config file path
    CONFIG_FILE="$1"
    shift
    CLI_MODE=false
else
    echo "ERROR: Unknown option '$1'"
    usage
    exit 1
fi

# ============================================================
# Environment setup
# ============================================================
echo "=========================================="
echo "PLIP Trajectory Analysis"
echo "=========================================="
if $CLI_MODE; then
    echo "Mode:       CLI (auto-detect chains)"
else
    echo "Mode:       Config file → $CONFIG_FILE"
fi
echo "Started at: $(date)"
echo ""

echo "Loading modules..."
module purge
module load python/3.12.3/gcc13.2/miniconda_24.4.0
module load pdb-tools/2.5.1
module load PLIP/3.12
module load amber/A24       # provides cpptraj

# Activate conda environment
CONDA_ENV="plip_analysis"
if conda env list | grep -q "^${CONDA_ENV} "; then
    echo "Activating conda environment: ${CONDA_ENV}"
    conda activate "${CONDA_ENV}"
else
    echo "WARNING: conda environment '${CONDA_ENV}' not found."
    echo "Creating it now..."
    conda create -n "${CONDA_ENV}" python=3.12 -y
    conda activate "${CONDA_ENV}"
    pip install pyyaml pandas matplotlib seaborn numpy tqdm
fi

echo ""
echo "Checking dependencies..."
which cpptraj  || { echo "ERROR: cpptraj not found!";  exit 1; }
which plip     || { echo "ERROR: plip not found!";     exit 1; }
which pdb_reres || { echo "ERROR: pdb-tools not found!"; exit 1; }
echo ""

# ============================================================
# Run analysis
# ============================================================
echo "Starting PLIP trajectory analysis..."
echo ""

if $CLI_MODE; then
    # Pass all remaining arguments directly to the Python script
    python plip_trajectory_analysis.py "$@"
else
    # Config file mode
    if [[ ! -f "$CONFIG_FILE" ]]; then
        echo "ERROR: Config file not found: $CONFIG_FILE"
        exit 1
    fi
    python plip_trajectory_analysis.py "$CONFIG_FILE"
fi

EXIT_CODE=$?

# ============================================================
# Done
# ============================================================
echo ""
echo "=========================================="
if [[ $EXIT_CODE -eq 0 ]]; then
    echo "Analysis completed successfully!"
else
    echo "Analysis FAILED (exit code: $EXIT_CODE)"
fi
echo "Finished at: $(date)"
echo "=========================================="

conda deactivate
module purge

exit $EXIT_CODE
