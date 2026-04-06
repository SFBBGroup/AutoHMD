#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to fix TER records in pdb4amber output files.
Removes redundant TER cards within the same chain while preserving
TER cards that separate different chains or mark true chain termini.
"""
import glob

def ter_fix():
    """
    Reads HIS_CYX_pdb4amber.pdb and removes redundant TER records.
    
    A TER record is considered redundant if:
    - The next line is an ATOM/HETATM record
    - AND it belongs to the same chain
    
    This preserves TER records that:
    - Separate different chains
    - Mark the end of a chain (before END or another chain)
    - Are at the end of the file
    """
    for ter in glob.glob("HIS_CYX_pdb4amber.pdb"):
        with open(ter, "r") as f_in, open("out_pdb4amber-addH.pdb", "w") as f_out:
            lines = f_in.readlines()
            corrected = []
            
            for i, line in enumerate(lines):
                # If it's not a TER record, keep it as-is:
                if not line.startswith("TER"):
                    corrected.append(line)
                    continue
                
                # It is TER - look ahead to the next line if there is one:
                if i + 1 < len(lines):
                    nxt = lines[i+1]
                else:
                    # No next line (EOF), so definitely keep this TER:
                    corrected.append(line)
                    continue
                
                # Only compare chain IDs if both lines are long enough:
                if len(line) > 21 and len(nxt) > 21:
                    chain_here = line[21]
                    chain_next = nxt[21]
                    
                    # If the next record is ATOM/HETATM in the same chain - skip this TER
                    if (nxt.startswith("ATOM") or nxt.startswith("HETATM")) and chain_here == chain_next:
                        continue
                
                # All other cases: keep the TER
                corrected.append(line)
            
            f_out.writelines(corrected)
            print(f"Processed: {ter}")
            print(f"Output written to: out_pdb4amber-addH.pdb")

if __name__ == "__main__":
    ter_fix()
