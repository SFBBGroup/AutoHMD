#!/usr/bin/env python3
"""
Script to detect gaps in PDB structures and insert TER cards.
A gap is defined when there is a difference > 1 in residue number
within the same chain.
IMPORTANT: Run BEFORE renumbering (pdb_reres) to detect original gaps.
"""
import sys
import re

def parse_pdb_line(line):
    """Extracts relevant information from an ATOM/HETATM line."""
    if not (line.startswith('ATOM') or line.startswith('HETATM')):
        return None
    
    try:
        record_type = line[0:6].strip()
        chain_id = line[21:22].strip()
        res_num = int(line[22:26].strip())
        res_name = line[17:20].strip()
        atom_num = int(line[6:11].strip())
        
        return {
            'record': record_type,
            'chain': chain_id,
            'resnum': res_num,
            'resname': res_name,
            'atomnum': atom_num
        }
    except (ValueError, IndexError):
        return None

def fix_gaps_with_ter(input_file, output_file, gap_threshold=1):
    """
    Reads a PDB file and inserts TER cards where there are gaps in numbering.
    
    Args:
        input_file: input PDB file
        output_file: output PDB file with inserted TERs
        gap_threshold: minimum difference to consider a gap (default: 1)
    """
    
    with open(input_file, 'r') as f:
        lines = f.readlines()
    
    output_lines = []
    prev_chain = None
    prev_resnum = None
    prev_atomnum = None
    gaps_found = []
    
    for line in lines:
        parsed = parse_pdb_line(line)
        
        if parsed:
            current_chain = parsed['chain']
            current_resnum = parsed['resnum']
            current_atomnum = parsed['atomnum']
            
            # Check if there is a gap
            if (prev_chain is not None and 
                prev_resnum is not None and 
                current_chain == prev_chain and
                current_resnum - prev_resnum > gap_threshold):
                
                # Insert TER with incremented atom number
                ter_atomnum = prev_atomnum + 1 if prev_atomnum else current_atomnum
                ter_line = f"TER   {ter_atomnum:5d}      {parsed['resname']} {current_chain}{prev_resnum:4d}\n"
                output_lines.append(ter_line)
                
                gap_info = f"Gap detected in chain {current_chain}: residue {prev_resnum} -> {current_resnum}"
                gaps_found.append(gap_info)
                print(gap_info)
            
            prev_chain = current_chain
            prev_resnum = current_resnum
            prev_atomnum = current_atomnum
        
        output_lines.append(line)
    
    # Write output file
    with open(output_file, 'w') as f:
        f.writelines(output_lines)
    
    print(f"\nFile processed: {output_file}")
    print(f"Total gaps found: {len(gaps_found)}")
    
    return len(gaps_found)

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python3 fix_gaps_with_ter.py <file.pdb> [output.pdb] [threshold]")
        print("\nExample:")
        print("  python3 fix_gaps_with_ter.py processed.pdb processed_fixed.pdb")
        print("  python3 fix_gaps_with_ter.py processed.pdb processed_fixed.pdb 2")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2] if len(sys.argv) > 2 else input_file.replace('.pdb', '_fixed.pdb')
    gap_threshold = int(sys.argv[3]) if len(sys.argv) > 3 else 1
    
    print(f"Processing: {input_file}")
    print(f"Threshold for gap: >{gap_threshold}")
    print("-" * 50)
    
    fix_gaps_with_ter(input_file, output_file, gap_threshold)
