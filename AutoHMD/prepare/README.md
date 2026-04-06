
# prepare/ — System Preparation

This folder contains the scripts and configuration file needed to convert a raw PDB file into a fully prepared, solvated, and ionised AMBER topology ready for Heated Molecular Dynamics simulation.

---


## Files

| File | Description |
|---|---|
| `config_hMD.conf` | User configuration file — edit before running |
| `prepare_autohMD.sh` | Master preparation script |
| `scripts/` | Auxiliary Python scripts called automatically (see `scripts/README.md`) |

---

## Configuration — `config_hMD.conf`

Edit this file before running `prepare_autohMD.sh`. It is the only file that needs to be modified for a standard run.

```bash
# Input PDB file name (without .pdb extension)
INPUT_PDB="your_protein"

# Protein force field (examples: ff19SB, ff14SB, fb15)
FORCEFIELD="ff19SB"

# Water model (examples: opc, tip3p, tip4pew, spce)
WATER_MODEL="opc"

# pH for protonation state calculation
PH="7.4"

# Solvation distance in Angstroms
SOLVATION_DISTANCE="12.0"

# Ion concentration in mol/L (physiological default: 0.15)
ION_CONCENTRATION="0.15"
```

### Parameter reference

| Parameter | Default | Description |
|---|---|---|
| `INPUT_PDB` | *(required)* | PDB filename without the `.pdb` extension. The file must be present in the `prepare/` directory. |
| `FORCEFIELD` | `ff19SB` | AMBER protein force field passed to `tleap`. |
| `WATER_MODEL` | `opc` | Explicit water model. Supported values: `opc`, `tip3p`, `tip4pew`, `spce`. |
| `PH` | `7.4` | pH value used by PDB2PQR/PROPKA to assign protonation states, histidine tautomers, and disulfide bridges. |
| `SOLVATION_DISTANCE` | `12.0` | Minimum distance in Å from any solute atom to the periodic box edge. |
| `ION_CONCENTRATION` | `0.15` | Target NaCl concentration in mol/L used to calculate the number of Na⁺ and Cl⁻ ions. |

---

## Usage

1. Place your PDB file in the `prepare/` directory.
2. Edit `config_hMD.conf` with at minimum the `INPUT_PDB` value.
3. Run the preparation script:

```bash
bash prepare_autohMD.sh
```

The script will stop with an informative error message if any step fails. All intermediate files are written to the current directory so that individual steps can be inspected or re-run if needed.

On successful completion, the three files required for simulation (`HMR_solvated_system.prmtop`, `solvated_system.rst7`, `solvated_system.pdb`) are automatically copied to `../submit/`.

---

## What `prepare_autohMD.sh` does

The script executes the following steps in order:

### Step 1 — Validation
Checks that `config_hMD.conf` exists, that `INPUT_PDB` has been set to a non-default value, and that the corresponding `.pdb` file is present. Checks that the `scripts/` directory exists with the required auxiliary scripts.

### Step 2 — PDB cleaning
Removes lines shorter than 26 characters (malformed or incomplete records) that would cause downstream tools to fail.

```
valid_fixed.pdb
```

### Step 3 — Gap detection and TER insertion
Calls `scripts/fix_gaps_with_ter.py` to detect residue numbering gaps within each chain and insert `TER` records at each gap. This must happen before renumbering so that the original residue numbers are used for gap detection.

```
valid_with_ter.pdb
```

### Step 4 — Residue renumbering
Uses `pdb_reres -1` and `pdb_fixinsert` (from PDBTools) to renumber all residues starting from 1 and resolve insertion codes.

```
renumbered.pdb
```

### Step 5 — Protonation state assignment
Runs PDB2PQR with the PROPKA titration method at the configured pH. Assigns histidine tautomers (HID/HIE/HIP), cysteine states (CYS/CYX), and other pH-dependent protonation states using the AMBER naming convention.

```
HIS_CYX_AMBER  (PQR output)
HIS_CYX.pdb    (PDB output)
```

### Step 6 — pdb4amber standardisation
Runs `pdb4amber` to standardise atom names, remove alternate conformations, and ensure AMBER compatibility.

```
HIS_CYX_pdb4amber.pdb
```

### Step 7 — TER record fix
Calls `scripts/pdb4amber_ter_fix.py` to remove spurious `TER` records that `pdb4amber` sometimes inserts within a single chain, which would cause `tleap` to misinterpret the topology.

```
out_pdb4amber-addH.pdb
```

### Step 8 — Hydrogen removal
Uses `reduce -Trim` to strip all hydrogens from the structure. This ensures that hydrogens are added consistently and correctly by `tleap` in a later step.

```
leap_input.pdb
```

### Step 9 — Terminal capping
Calls `scripts/capping.py` to apply ACE (acetyl) caps at N-termini and NME (N-methyl amide) caps at C-termini of every chain fragment defined by `TER` records. This prevents charged termini from introducing artefacts during simulation.

```
leap_input.pdb  (updated in place)
```

### Step 10 — Initial solvation
Generates `tleap_initial.in` using the configured force field, water model, and solvation distance, then runs `tleap` to produce an initial solvated system. This run is used to extract the system volume and net charge for ion calculation.

```
tleap_initial.in
tleap.log
solvated_system.prmtop  (initial, before ion correction)
solvated_system.rst7
nowat_system.prmtop / nowat_system.rst7
```

### Step 11 — Ion concentration calculation
Calls `scripts/ion_concentration.py`, which reads the box volume and net charge from `tleap.log` and calculates the number of Na⁺ and Cl⁻ ions required to reach the target salt concentration while neutralising the system. Writes the final `tleap.in`.

```
tleap.in
```

### Step 12 — Final solvation
Runs `tleap` with the corrected ion counts to produce the final solvated and ionised system.

```
tleap_final.log
solvated_system.prmtop  (final)
solvated_system.rst7    (final)
solvated_system.pdb
```

### Step 13 — Hydrogen Mass Repartitioning (HMR)
Runs `parmed` to repartition hydrogen masses, transferring mass from heavy atoms to bonded hydrogens. HMR allows the use of a 4 fs timestep in production without loss of accuracy.

```
HMR_solvated_system.prmtop
parmed.log
```

---

