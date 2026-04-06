#!/usr/bin/env python3
"""
Script to detect gaps in PDB structures and insert TER records to break chains and prevent artificial bonds in molecular dynamics simulations.

Usage:
    python detect_and_fix_gaps.py input.pdb --check-only          # Check
    python detect_and_fix_gaps.py input.pdb -o output.pdb         # Fix
"""

import sys
import argparse
import math
from collections import defaultdict

def parse_pdb_atoms(filename):
    """Reads the PDB file and extracts atom information."""
    atoms = []
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith(('ATOM', 'HETATM')):
                try:
                    atom_info = {
                        'line': line,
                        'chain': line[21:22].strip() or 'A',
                        'resnum': int(line[22:26].strip()),
                        'resname': line[17:20].strip(),
                        'atom_name': line[12:16].strip(),
                        'x': float(line[30:38].strip()),
                        'y': float(line[38:46].strip()),
                        'z': float(line[46:54].strip())
                    }
                    atoms.append(atom_info)
                except (ValueError, IndexError):
                    continue
    return atoms

def calculate_distance(atom1, atom2):
    """Calculates the Euclidean distance between two atoms."""
    dx = atom1['x'] - atom2['x']
    dy = atom1['y'] - atom2['y']
    dz = atom1['z'] - atom2['z']
    return math.sqrt(dx*dx + dy*dy + dz*dz)

def detect_gaps(atoms, distance_threshold=4.0):
    """
    Detects gaps by analyzing C–N distances between consecutive residues.
    
    A normal peptide bond has a length of approximately ~1.33 Å. Structural gaps present much larger C–N distances (typically > 10 Å).
    
    Args:
        atoms: List of atoms from the PDB file
        distance_threshold: Maximum C–N distance for a valid bond (Å)
    
    Returns:
        List of dictionaries containing information about the detected gaps
    """
    gaps = []
    chains = defaultdict(list)
    
    # Organize atoms by chain
    for atom in atoms:
        chains[atom['chain']].append(atom)
    
    # Analyze each chain separately
    for chain_id, chain_atoms in chains.items():
        # Sort by residue number
        chain_atoms.sort(key=lambda x: x['resnum'])
        
        # Group atoms by residue
        residues = {}
        for atom in chain_atoms:
            resnum = atom['resnum']
            if resnum not in residues:
                residues[resnum] = {'atoms': [], 'resname': atom['resname']}
            residues[resnum]['atoms'].append(atom)
        
        res_numbers = sorted(residues.keys())
        
        # Check gaps between consecutive residues in numbering
        for i in range(len(res_numbers) - 1):
            curr_res = res_numbers[i]
            next_res = res_numbers[i + 1]
            
            # Find C atom (carbonyl) of current residue
            c_atom = None
            for atom in residues[curr_res]['atoms']:
                if atom['atom_name'] == 'C':
                    c_atom = atom
                    break
            
            # Find N atom (amide) of next residue
            n_atom = None
            for atom in residues[next_res]['atoms']:
                if atom['atom_name'] == 'N':
                    n_atom = atom
                    break
            
            # If both atoms found, check distance
            if c_atom and n_atom:
                distance = calculate_distance(c_atom, n_atom)
                sequence_gap = next_res - curr_res
                
                # Gap detected: large distance OR jump in numbering
                if distance > distance_threshold or sequence_gap > 1:
                    gap_info = {
                        'chain': chain_id,
                        'res_before': curr_res,
                        'res_after': next_res,
                        'resname_before': residues[curr_res]['resname'],
                        'resname_after': residues[next_res]['resname'],
                        'distance': distance,
                        'sequence_gap': sequence_gap,
                        'missing_residues': sequence_gap - 1
                    }
                    gaps.append(gap_info)
    
    return gaps

def print_gap_report(gaps):
    """Prints detailed report of found gaps."""
    if not gaps:
        print("✓ No gaps detected in structure!")
        return
    
    print(f"\n{'='*75}")
    print(f"GAPS DETECTED: {len(gaps)}")
    print(f"{'='*75}\n")
    
    for i, gap in enumerate(gaps, 1):
        print(f"Gap #{i}:")
        print(f"  Chain: {gap['chain']}")
        print(f"  Position: {gap['resname_before']}{gap['res_before']} → "
              f"{gap['resname_after']}{gap['res_after']}")
        print(f"  C-N Distance: {gap['distance']:.2f} Å", end='')
        
        if gap['distance'] > 10:
            print(" [STRUCTURAL GAP]")
        elif gap['distance'] > 4:
            print(" [Suspicious distance]")
        else:
            print()
        
        print(f"  Missing residues: {gap['missing_residues']}")
        print()

def insert_ter_cards(input_file, output_file, gaps):
    """
    Inserts TER cards after residues before gaps.
    
    This informs tleap that chains should be treated
    separately, avoiding artificial bonds.
    """
    # Create set of positions where to insert TER
    ter_positions = {}
    for gap in gaps:
        key = (gap['chain'], gap['res_before'])
        ter_positions[key] = gap
    
    lines_written = 0
    ter_count = 0
    
    with open(input_file, 'r') as fin, open(output_file, 'w') as fout:
        prev_line_was_atom = False
        prev_chain = None
        prev_resnum = None
        ter_inserted = set()
        
        for line in fin:
            if line.startswith(('ATOM', 'HETATM')):
                chain = line[21:22].strip() or 'A'
                resnum = int(line[22:26].strip())
                
                # Check if TER needs to be inserted
                key = (chain, prev_resnum)
                if (key in ter_positions and 
                    key not in ter_inserted and 
                    resnum != prev_resnum):  # Changed residue
                    
                    fout.write("TER\n")
                    ter_inserted.add(key)
                    ter_count += 1
                    gap = ter_positions[key]
                    
                    print(f"✓ TER inserted after: {gap['resname_before']}{gap['res_before']} "
                          f"(chain {chain})")
                    print(f"  Next segment starts at: {gap['resname_after']}{gap['res_after']}")
                
                fout.write(line)
                lines_written += 1
                prev_line_was_atom = True
                prev_chain = chain
                prev_resnum = resnum
            
            elif line.startswith('TER'):
                # Keep existing TER records
                if prev_line_was_atom:
                    fout.write(line)
                    prev_line_was_atom = False
            
            else:
                # Other lines (HEADER, REMARK, CONECT, END, etc)
                fout.write(line)
                prev_line_was_atom = False
        
        # Add final TER if ended with ATOM
        if prev_line_was_atom:
            fout.write("TER\n")
    
    return ter_count, lines_written

def main():
    parser = argparse.ArgumentParser(
        description='Detects structural gaps and inserts TER cards for MD',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Only check for gaps
  python detect_and_fix_gaps.py input.pdb --check-only
  
  # Fix and save
  python detect_and_fix_gaps.py input.pdb -o output.pdb
  
  # Adjust distance threshold
  python detect_and_fix_gaps.py input.pdb -o output.pdb --distance-threshold 3.5

Notes:
  - Normal peptide bond: ~1.33 Å
  - Default threshold: 4.0 Å (conservative)
  - Typical structural gaps: > 10 Å
        """
    )
    
    parser.add_argument('input', help='Input PDB file')
    parser.add_argument('-o', '--output', 
                        help='Output PDB file (default: input_fixed.pdb)')
    parser.add_argument('--distance-threshold', type=float, default=4.0,
                        help='Maximum C-N distance for valid bond in Å (default: 4.0)')
    parser.add_argument('--check-only', action='store_true',
                        help='Only detect gaps without modifying file')
    parser.add_argument('-q', '--quiet', action='store_true',
                        help='Quiet mode (errors only)')
    
    args = parser.parse_args()
    
    # Define output filename if not provided
    if not args.output:
        base = args.input.rsplit('.', 1)[0]
        args.output = f"{base}_fixed.pdb"
    
    if not args.quiet:
        print(f"\n{'='*75}")
        print(f"GAP ANALYSIS - {args.input}")
        print(f"{'='*75}")
    
    # Parse PDB
    try:
        atoms = parse_pdb_atoms(args.input)
    except FileNotFoundError:
        print(f"\n✗ Error: File '{args.input}' not found!")
        return 1
    
    if not atoms:
        print("\n✗ Error: No atoms found in file!")
        return 1
    
    if not args.quiet:
        chains = set(a['chain'] for a in atoms)
        print(f"\nAtoms read: {len(atoms)}")
        print(f"Chains: {', '.join(sorted(chains))}")
        print(f"Distance threshold: {args.distance_threshold} Å")
    
    # Detect gaps
    gaps = detect_gaps(atoms, distance_threshold=args.distance_threshold)
    
    # Print report
    if not args.quiet:
        print_gap_report(gaps)
    
    # If no gaps, just copy file
    if not gaps:
        if not args.check_only:
            import shutil
            shutil.copy(args.input, args.output)
            if not args.quiet:
                print(f"\n✓ File copied to: {args.output}")
                print("  (no modifications needed)\n")
        return 0
    
    # Check-only mode: don't modify file
    if args.check_only:
        print("\n[--check-only mode: no corrections applied]")
        print(f"To fix, run:")
        print(f"  python {sys.argv[0]} {args.input} -o {args.output}\n")
        return 0
    
    # Insert TER cards
    if not args.quiet:
        print(f"\n{'='*75}")
        print("FIXING GAPS")
        print(f"{'='*75}\n")
    
    ter_count, lines_written = insert_ter_cards(args.input, args.output, gaps)
    
    if not args.quiet:
        print(f"\n{'='*75}")
        print("✓ CORRECTION COMPLETED!")
        print(f"{'='*75}")
        print(f"\nOutput file: {args.output}")
        print(f"Lines written: {lines_written}")
        print(f"TER cards inserted: {ter_count}")
        print(f"\nNow tleap will treat each segment as an independent chain.")
        print("No more long bond warnings!\n")
    
    return 0

if __name__ == '__main__':
    sys.exit(main())
