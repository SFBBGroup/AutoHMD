#!/bin/bash
# Automated script for protein system preparation for Heated Molecular Dynamics

set -e  # Stop execution on error

# Determine script directory and scripts folder
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPTS_DIR="${SCRIPT_DIR}/scripts"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Logging functions
log_info() {
    echo -e "${GREEN}[INFO]${NC} $1"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

log_warn() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

# Check if scripts directory exists
if [ ! -d "$SCRIPTS_DIR" ]; then
    log_error "Scripts directory not found: $SCRIPTS_DIR"
    log_info "Please create a 'scripts' folder with the required Python scripts."
    exit 1
fi

# Check if configuration file exists
if [ ! -f "config_hMD.conf" ]; then
    log_error "config_hMD.conf file not found!"
    log_info "Creating template configuration file..."
    cat > config_hMD.conf << 'EOF'
# Configuration for Heated Molecular Dynamics Preparation
# Edit the values below as needed

# Input PDB file name (without .pdb extension)
INPUT_PDB="my_protein"

# Protein force field (examples: ff19SB, ff14SB, fb15)
FORCEFIELD="ff19SB"

# Water model (examples: opc, tip3p, tip4pew, spce)
WATER_MODEL="opc"

# pH for protonation state calculation
PH="7.4"

# Solvation distance in Angstroms (default: 12.0)
SOLVATION_DISTANCE="12.0"
EOF
    log_info "config_hMD.conf file created. Please edit it and run the script again."
    exit 1
fi

# Load configurations
log_info "Loading configurations..."
source config_hMD.conf

# Validate configurations
if [ -z "$INPUT_PDB" ] || [ "$INPUT_PDB" == "my_protein" ]; then
    log_error "INPUT_PDB was not configured in config_hMD.conf"
    exit 1
fi

if [ ! -f "${INPUT_PDB}.pdb" ]; then
    log_error "File ${INPUT_PDB}.pdb not found!"
    exit 1
fi

log_info "Configurations loaded:"
log_info "  - Input PDB: ${INPUT_PDB}.pdb"
log_info "  - Force field: ${FORCEFIELD}"
log_info "  - Water model: ${WATER_MODEL}"
log_info "  - pH: ${PH}"
log_info "  - Solvation distance: ${SOLVATION_DISTANCE} Å"
log_info "  - Scripts directory: ${SCRIPTS_DIR}"

# Load modules
log_info "Loading required modules..."
module purge
module load amber/A24 pdb2pqr/3.6.1 pdb-tools/2.5.1

# Process PDB
log_info "Processing PDB file..."
file_pdb="${INPUT_PDB}.pdb"

# Remove malformed lines
log_info "Removing malformed lines..."
awk 'length($0) >= 26' "$file_pdb" > valid_fixed.pdb

# Insert TER records
log_info "Fixing gaps with TER records..."
if [ -f "${SCRIPTS_DIR}/fix_gaps_with_ter.py" ]; then
    python3 "${SCRIPTS_DIR}/fix_gaps_with_ter.py" valid_fixed.pdb valid_with_ter.pdb
else
    log_warn "fix_gaps_with_ter.py not found in ${SCRIPTS_DIR}, skipping this step..."
    cp valid_fixed.pdb valid_with_ter.pdb
fi

# Renumber residues
log_info "Renumbering residues..."
pdb_reres -1 valid_with_ter.pdb | pdb_fixinsert > renumbered.pdb

# Calculate protonation states
log_info "Calculating protonation states at pH ${PH}..."
pdb2pqr30 --ff='AMBER' renumbered.pdb \
    --titration-state-method=propka \
    --with-ph "${PH}" \
    --ffout='AMBER' HIS_CYX_AMBER \
    --pdb-output HIS_CYX.pdb

# Process with pdb4amber
log_info "Standardizing structure with pdb4amber..."
pdb4amber HIS_CYX.pdb -o HIS_CYX_pdb4amber.pdb

# Fix TER records after pdb4amber
log_info "Fixing TER records..."
if [ -f "${SCRIPTS_DIR}/pdb4amber_ter_fix.py" ]; then
    python3 "${SCRIPTS_DIR}/pdb4amber_ter_fix.py"
else
    log_warn "pdb4amber_ter_fix.py not found in ${SCRIPTS_DIR}, skipping this step..."
    cp HIS_CYX_pdb4amber.pdb out_pdb4amber-addH.pdb
fi

# Remove hydrogens
log_info "Removing hydrogens for consistent re-addition..."
reduce -Trim out_pdb4amber-addH.pdb > leap_input.pdb 2>&1 || true
sleep 1  # Give time for file to be written
if [ ! -s "leap_input.pdb" ]; then
    log_error "reduce failed to generate leap_input.pdb!"
    exit 1
fi
log_info "Hydrogens removed successfully"

# Cap protein termini
log_info "Capping protein C- and N-terminals..."
if [ -f "${SCRIPTS_DIR}/capping.py" ]; then
    python3 "${SCRIPTS_DIR}/capping.py" leap_input.pdb || true
    
    if [ -f "TEMP_capped.pdb" ]; then
        mv TEMP_capped.pdb leap_input.pdb
        log_info "Protein termini capped successfully!"
    else
        log_warn "TEMP_capped.pdb not generated, continuing without capping..."
    fi
else
    log_warn "capping.py not found in ${SCRIPTS_DIR}, skipping terminal capping..."
    log_warn "Protein will have charged termini (NH3+ and COO-)"
fi

# Verify leap_input.pdb exists and is not empty
log_info "Verifying leap_input.pdb..."
if [ ! -s "leap_input.pdb" ]; then
    log_error "leap_input.pdb is missing or empty!"
    exit 1
fi

leap_lines=$(wc -l < leap_input.pdb)
log_info "leap_input.pdb verified ($leap_lines lines)"
log_info "Proceeding to tleap file generation..."

# Determine water model leaprc
case "${WATER_MODEL,,}" in
    opc)
        water_leaprc="leaprc.water.opc"
        water_box="OPCBOX"
        ;;
    tip3p)
        water_leaprc="leaprc.water.tip3p"
        water_box="TIP3PBOX"
        ;;
    tip4pew)
        water_leaprc="leaprc.water.tip4pew"
        water_box="TIP4PEWBOX"
        ;;
    spce)
        water_leaprc="leaprc.water.spce"
        water_box="SPCBOX"
        ;;
    *)
        log_error "Water model '${WATER_MODEL}' not recognized!"
        log_info "Supported models: opc, tip3p, tip4pew, spce"
        exit 1
        ;;
esac

# Create tleap_initial.in file
log_info "Generating tleap_initial.in file..."
cat > tleap_initial.in << EOF
# tleap file automatically generated by step1_autohMD.sh
# Settings: FF=${FORCEFIELD}, Water=${WATER_MODEL}, pH=${PH}

source leaprc.protein.${FORCEFIELD}
source ${water_leaprc}

pd = loadpdb leap_input.pdb
check pd

# Save system without water
saveamberparm pd nowat_system.prmtop nowat_system.rst7

# Check charge and neutralize system
charge pd
addions pd Na+ 0
addions pd Cl- 0

# Solvate with octahedral box
solvateOct pd ${water_box} ${SOLVATION_DISTANCE}

# Save solvated system
saveamberparm pd solvated_system.prmtop solvated_system.rst7
savePdb pd solvated_system.pdb

quit
EOF

log_info "tleap_initial.in file created successfully!"
log_info ""
log_info "=========================================="
log_info "Step 1: Initial preparation completed!"
log_info "=========================================="
log_info ""
log_info "Running initial tleap to get system volume and charge..."
tleap -f tleap_initial.in > tleap.log 2>&1

if [ ! -f "solvated_system.prmtop" ]; then
    log_error "Initial tleap failed! Check tleap.log for errors."
    exit 1
fi

# Extract and display 'check pd' output from tleap.log
log_info ""
log_info "=========================================="
log_info "System Check Results (from tleap):"
log_info "=========================================="
if grep -q "Checking 'pd'" tleap.log; then
    # Extract from "Checking 'pd'" to "Unit is OK" or end of check section
    awk '/Checking .pd/{flag=1} flag{print} /Unit is OK|Unit is not OK/{flag=0}' tleap.log | while read -r line; do
        if [[ "$line" == *"Warning"* ]]; then
            echo -e "${YELLOW}${line}${NC}"
        elif [[ "$line" == *"Unit is OK"* ]]; then
            echo -e "${GREEN}${line}${NC}"
        elif [[ "$line" == *"Unit is not OK"* ]]; then
            echo -e "${RED}${line}${NC}"
        else
            echo "$line"
        fi
    done
else
    log_warn "Could not find 'check pd' output in tleap.log"
fi
log_info "=========================================="
log_info ""

log_info "Initial solvation completed successfully!"

log_info ""
log_info "=========================================="
log_info "Step 2: Calculating ion concentration..."
log_info "=========================================="
python3 "${SCRIPTS_DIR}/ion_concentration.py"

if [ ! -f "tleap.in" ]; then
    log_error "ion_concentration.py failed to generate tleap.in!"
    exit 1
fi

log_info ""
log_info "Running final tleap with calculated ion concentrations..."
tleap -f tleap.in > tleap_final.log 2>&1

if [ ! -f "solvated_system.prmtop" ] || [ ! -f "solvated_system.rst7" ]; then
    log_error "Final tleap failed! Check tleap_final.log for errors."
    exit 1
fi
log_info "System solvated and ionized successfully!"

log_info ""
log_info "=========================================="
log_info "Step 3: Applying Hydrogen Mass Repartitioning..."
log_info "=========================================="

# Create parmed.in file for HMR
cat > parmed.in << 'EOF'
parm solvated_system.prmtop
HMassRepartition
outparm HMR_solvated_system.prmtop
quit
EOF

log_info "Running ParmEd for HMR..."
parmed -i parmed.in > parmed.log 2>&1

if [ ! -f "HMR_solvated_system.prmtop" ]; then
    log_error "ParmEd HMR failed! Check parmed.log for errors."
    exit 1
fi
log_info "Hydrogen Mass Repartitioning completed successfully!"

cp HMR_solvated_system.prmtop solvated_system.rst7 solvated_system.pdb ../submit

log_info ""
log_info "=========================================="
log_info "ALL STEPS COMPLETED SUCCESSFULLY!"
log_info "=========================================="
log_info ""
log_info "Final output files transfered to submission folder:"
log_info "  - HMR_solvated_system.prmtop (topology with HMR)"
log_info "  - solvated_system.rst7 (coordinates)"
log_info "  - solvated_system.pdb (for visualization)"
log_info ""
log_info "System is ready for Molecular Dynamics!"
log_info ""
log_info "Log files created:"
log_info "  - tleap.log (initial solvation)"
log_info "  - tleap_final.log (final system with ions)"
log_info "  - parmed.log (HMR application)"
