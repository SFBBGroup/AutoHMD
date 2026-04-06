#!/usr/bin/env python3
"""
Calculate ion concentration and generate tleap input file for MD system preparation.
"""
import re
import sys

def read_config(config_file="config_hMD.conf"):
    """Read configuration from config file."""
    config = {}
    try:
        with open(config_file) as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith('#'):
                    if '=' in line:
                        key, value = line.split('=', 1)
                        config[key.strip()] = value.strip().strip('"').strip("'")
    except FileNotFoundError:
        print(f"Error: {config_file} not found!")
        sys.exit(1)
    return config

# Read configuration
config = read_config()
forcefield = config.get('FORCEFIELD', 'ff19SB')
water_model = config.get('WATER_MODEL', 'opc').lower()
solvation_distance = config.get('SOLVATION_DISTANCE', '12.0')
ion_concentration = float(config.get('ION_CONCENTRATION', '0.15'))

# Determine water model parameters
water_leaprc_map = {
    'opc': ('leaprc.water.opc', 'OPCBOX'),
    'tip3p': ('leaprc.water.tip3p', 'TIP3PBOX'),
    'tip4pew': ('leaprc.water.tip4pew', 'TIP4PEWBOX'),
    'spce': ('leaprc.water.spce', 'SPCBOX')
}

if water_model not in water_leaprc_map:
    print(f"Error: Unsupported water model '{water_model}'")
    sys.exit(1)

water_leaprc, water_box = water_leaprc_map[water_model]

# Extract values from tleap.log
try:
    with open("tleap.log") as f:
        content = f.read()
        # Capture charge value (including 0.000000 or negative values)
        charge_match = re.search(r"Total unperturbed charge:\s*(-?\d+\.\d+)", content)
        # Capture volume
        volume_match = re.search(r"Volume:\s*([0-9.]+)", content)
        
        if charge_match:
            charge_protein = float(charge_match.group(1))
        else:
            charge_protein = 0.0  # safe default if not found
        
        if volume_match:
            volume_angstrom3 = float(volume_match.group(1))
        else:
            raise ValueError("Volume not found in tleap.log file")
except FileNotFoundError:
    print("Error: tleap.log not found!")
    print("Please run: tleap -f tleap_initial.in > tleap.log 2>&1")
    sys.exit(1)

# Parameters
avogadro_number = 6.022e23
water_molecule_volume_A3 = 30

# Calculations
num_water_molecules = volume_angstrom3 / water_molecule_volume_A3
num_ions_expected = num_water_molecules * ion_concentration / 56
num_na = round(num_ions_expected - charge_protein / 2)
num_cl = round(num_ions_expected + charge_protein / 2)
mol_na = num_na / avogadro_number
mol_cl = num_cl / avogadro_number
volume_L = volume_angstrom3 * 1e-27  # converting Å³ to L
molarity_na = mol_na / volume_L
molarity_cl = mol_cl / volume_L
molarity_total = (mol_na + mol_cl) / volume_L

# Results
print(f"Number of Na+ ions: {num_na}")
print(f"Number of Cl- ions: {num_cl}")
print(f"Molarity of Na+: {molarity_na:.3f} M")
print(f"Molarity of Cl-: {molarity_cl:.3f} M")
print(f"Total molarity: {molarity_total:.3f} M")
print()

# Generate tleap.in file
leap_content = f"""source leaprc.protein.{forcefield}
source {water_leaprc}
pd = loadpdb leap_input.pdb
check pd
saveamberparm pd nowat_system.prmtop nowat_system.rst7
charge pd
addions pd Na+ {num_na}
addions pd Cl- {num_cl}
solvateOct pd {water_box} {solvation_distance}
saveamberparm pd solvated_system.prmtop solvated_system.rst7
savePdb pd solvated_system.pdb
quit
"""

with open("tleap.in", "w") as f:
    f.write(leap_content)

print("File tleap.in successfully generated!")
print(f"\nContent of tleap.in:")
print(leap_content)
