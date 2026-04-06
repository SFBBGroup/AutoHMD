#!/usr/bin/env python3
"""
Script to cap N and C terminals with ACE and NME, preserving TER cards.
Caps ALL fragments created by TERs (gaps).
FIXED: Correct atom order in ACE (CH3, C, O)
"""

from Bio.PDB import PDBParser, PDBIO, Select
from Bio.PDB.Atom import Atom
import re

def detect_fragments_with_ter(pdb_file):
    """
    Detects fragments separated by TER cards.
    Returns: 
    - ter_positions: {(chain_id, res_num): True}
    - fragments: {chain_id: [(start_res, end_res), ...]}
    """
    ter_positions = {}
    fragments = {}
    
    current_chain = None
    fragment_start = None
    fragment_end = None
    prev_line_was_ter = False
    
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                try:
                    chain = line[21:22].strip()
                    resnum = int(line[22:26].strip())
                    
                    if current_chain != chain:
                        # New chain - save previous fragment if it exists
                        if current_chain and fragment_start is not None:
                            if current_chain not in fragments:
                                fragments[current_chain] = []
                            fragments[current_chain].append((fragment_start, fragment_end))
                        
                        # Start new chain
                        current_chain = chain
                        fragment_start = resnum
                        fragment_end = resnum
                        prev_line_was_ter = False
                    
                    elif prev_line_was_ter:
                        # After a TER, start new fragment in same chain
                        fragment_start = resnum
                        fragment_end = resnum
                        prev_line_was_ter = False
                    
                    else:
                        # Continuation of current fragment
                        fragment_end = resnum
                    
                except (ValueError, IndexError):
                    continue
            
            elif line.startswith('TER'):
                # TER found - mark end of current fragment
                if current_chain and fragment_start is not None:
                    if current_chain not in fragments:
                        fragments[current_chain] = []
                    fragments[current_chain].append((fragment_start, fragment_end))
                    ter_positions[(current_chain, fragment_end)] = True
                    print(f"TER detected: Chain {current_chain}, Residue {fragment_end}")
                
                # Mark that next ATOM starts new fragment
                prev_line_was_ter = True
                fragment_start = None
    
    # Save last fragment if it exists
    if current_chain and fragment_start is not None:
        if current_chain not in fragments:
            fragments[current_chain] = []
        fragments[current_chain].append((fragment_start, fragment_end))
    
    return ter_positions, fragments

def apply_capping(input_pdb, output_pdb, chains_ace_nme=None, chains_nme_only=None):
    """
    Applies ACE/NME capping to ALL fragments and reinserts TER cards.
    
    Args:
        input_pdb: input PDB file
        output_pdb: output PDB file
        chains_ace_nme: list of chains to cap N-term (ACE) and C-term (NME)
        chains_nme_only: list of chains to cap only C-term (NME)
    """
    
    # 1. Detect TERs and fragments
    ter_positions, fragments = detect_fragments_with_ter(input_pdb)
    
    print("\n=== FRAGMENTS DETECTED ===")
    for chain, frags in fragments.items():
        print(f"Chain {chain}: {len(frags)} fragment(s)")
        for i, (start, end) in enumerate(frags, 1):
            print(f"  Fragment {i}: residues {start}-{end}")
    
    # 2. Process structure with BioPython
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("model", input_pdb)
    
    # Define default chains if not specified
    all_chains = {chain.id for model in structure for chain in model}
    
    if chains_ace_nme is None:
        chains_ace_nme = [c for c in all_chains if c not in ("H", "L")]
    
    if chains_nme_only is None:
        chains_nme_only = ["H", "L"]
    
    print(f"\n=== CAPPING CONFIGURATION ===")
    print(f"Chains ACE+NME (all terminals): {chains_ace_nme}")
    print(f"Chains NME only (C-terminal): {chains_nme_only}")
    
    # 3. Apply capping to ALL fragments
    print(f"\n=== APPLYING CAPPING ===")
    
    for model in structure:
        for chain in model:
            chain_id = chain.id
            residues = list(chain.get_residues())
            
            if not residues:
                continue
            
            # Get fragments for this chain
            chain_fragments = fragments.get(chain_id, [])
            
            if chain_id in chains_ace_nme:
                # === Cap ALL terminals with ACE and NME ===
                for frag_idx, (start_res, end_res) in enumerate(chain_fragments):
                    # Find residues in fragment
                    frag_residues = [r for r in residues 
                                    if start_res <= r.id[1] <= end_res]
                    
                    if not frag_residues:
                        continue
                    
                    first = frag_residues[0]
                    last = frag_residues[-1]
                    
                    print(f"  Chain {chain_id}, Fragment {frag_idx+1} ({start_res}-{end_res}):")
                    
                    # === CREATE ACE AT N-TERMINAL WITH CORRECT ORDER ===
                    first.resname = "ACE"
                    
                    # Save original atom information
                    original_atoms = {}
                    for atom in list(first.get_atoms()):
                        name = atom.get_name()
                        if name in ("CA", "C", "O"):
                            original_atoms[name] = {
                                'coord': atom.coord.copy(),
                                'bfactor': atom.bfactor,
                                'occupancy': atom.occupancy,
                                'element': atom.element,
                                'serial': atom.serial_number
                            }
                    
                    # Remove ALL atoms from residue
                    for atom in list(first.get_atoms()):
                        first.detach_child(atom.id)
                    
                    # Recreate atoms in correct order: CH3, C, O
                    if "CA" in original_atoms:
                        info = original_atoms["CA"]
                        ch3 = Atom(
                            name="CH3",
                            coord=info['coord'],
                            bfactor=info['bfactor'],
                            occupancy=info['occupancy'],
                            altloc=' ',
                            fullname=" CH3",
                            serial_number=info['serial'],
                            element=info['element']
                        )
                        first.add(ch3)
                    
                    if "C" in original_atoms:
                        info = original_atoms["C"]
                        c = Atom(
                            name="C",
                            coord=info['coord'],
                            bfactor=info['bfactor'],
                            occupancy=info['occupancy'],
                            altloc=' ',
                            fullname=" C  ",
                            serial_number=info['serial'],
                            element=info['element']
                        )
                        first.add(c)
                    
                    if "O" in original_atoms:
                        info = original_atoms["O"]
                        o = Atom(
                            name="O",
                            coord=info['coord'],
                            bfactor=info['bfactor'],
                            occupancy=info['occupancy'],
                            altloc=' ',
                            fullname=" O  ",
                            serial_number=info['serial'],
                            element=info['element']
                        )
                        first.add(o)
                    
                    print(f"    N-term: ACE (res {first.id[1]}) - order: CH3, C, O")
                    
                    # === NME AT C-TERMINAL ===
                    last.resname = "NME"
                    atoms_to_remove = [a for a in last if a.get_name() != "N"]
                    for atom in atoms_to_remove:
                        last.detach_child(atom.id)
                    print(f"    C-term: NME (res {last.id[1]})")
            
            elif chain_id in chains_nme_only:
                # === Cap only C-terminals with NME ===
                for frag_idx, (start_res, end_res) in enumerate(chain_fragments):
                    frag_residues = [r for r in residues 
                                    if start_res <= r.id[1] <= end_res]
                    
                    if not frag_residues:
                        continue
                    
                    last = frag_residues[-1]
                    
                    print(f"  Chain {chain_id}, Fragment {frag_idx+1} ({start_res}-{end_res}):")
                    
                    # NME at C-terminal
                    last.resname = "NME"
                    atoms_to_remove = [a for a in last if a.get_name() != "N"]
                    for atom in atoms_to_remove:
                        last.detach_child(atom.id)
                    print(f"    C-term: NME (res {last.id[1]})")
    
    # 4. Save temporary structure
    temp_file = "temp_no_ter.pdb"
    io = PDBIO()
    io.set_structure(structure)
    io.save(temp_file)
    
    # 5. Reinsert TER cards at correct positions
    reinsert_ter_cards(temp_file, output_pdb, ter_positions)
    
    # Clean up temporary file
    import os
    os.remove(temp_file)
    
    print(f"\n✓ File saved: {output_pdb}")
    print(f"✓ TER cards preserved: {len(ter_positions)}")

def reinsert_ter_cards(input_pdb, output_pdb, ter_positions):
    """
    Reads the processed PDB and reinserts TER cards at original positions.
    REMOVES all TER and END lines from temporary file.
    """
    with open(input_pdb, 'r') as f:
        lines = f.readlines()
    
    output_lines = []
    prev_chain = None
    prev_resnum = None
    prev_resname = None
    atom_counter = 0
    
    for line in lines:
        # SKIP original TER and END lines - we'll insert our own TERs
        if line.startswith('TER') or line.startswith('END'):
            continue
        
        if line.startswith('ATOM') or line.startswith('HETATM'):
            try:
                chain = line[21:22].strip()
                resnum = int(line[22:26].strip())
                resname = line[17:20].strip()
                atom_counter = int(line[6:11].strip())
                
                # Check if TER should be inserted BEFORE this line
                if prev_chain and prev_resnum and (prev_chain, prev_resnum) in ter_positions:
                    ter_line = f"TER   {atom_counter:5d}      {prev_resname} {prev_chain}{prev_resnum:4d}\n"
                    output_lines.append(ter_line)
                
                prev_chain = chain
                prev_resnum = resnum
                prev_resname = resname
                
            except (ValueError, IndexError):
                pass
        
        output_lines.append(line)
    
    # Final TER if necessary
    if prev_chain and prev_resnum and (prev_chain, prev_resnum) in ter_positions:
        ter_line = f"TER   {atom_counter+1:5d}      {prev_resname} {prev_chain}{prev_resnum:4d}\n"
        output_lines.append(ter_line)
    
    # Add final END
    output_lines.append("END\n")
    
    # Write final file
    with open(output_pdb, 'w') as f:
        f.writelines(output_lines)

if __name__ == "__main__":
    import sys
    
    if len(sys.argv) < 2:
        print("Usage: python3 capping_with_ter.py <input.pdb> [output.pdb]")
        print("\nExample:")
        print("  python3 capping_with_ter.py processed.pdb TEMP_capped.pdb")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2] if len(sys.argv) > 2 else "TEMP_capped.pdb"
    
    print(f"Processing: {input_file}")
    print("=" * 60)
    
    apply_capping(input_file, output_file)
