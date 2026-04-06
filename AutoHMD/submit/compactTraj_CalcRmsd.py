#!/usr/bin/env python3
"""
Script to compact trajectories and automatically calculate iRMSD.
Executes two steps:
1. Compaction: removes water and ions from trajectories
2. Analysis: calculates iRMSD on the compacted trajectory

DETECTION LOGIC:
- Antigen: residues that start with ACE (can have multiple fragments)
- Antibody: residues that start with any other residue (always 2 chains with NME)
"""

import sys
import argparse
import subprocess
import os

def extract_residue_ranges(pdb_file):
    """
    Extracts residue ranges for antigen and antibody from a PDB file.
    
    DEFINITIVE LOGIC:
    - If a fragment starts with ACE -> it's ANTIGEN
    - If a fragment DOESN'T start with ACE -> it's ANTIBODY
    - Antibody always has 2 chains (first 2 fragments without ACE)
    - Antigen can have multiple fragments (all start with ACE)
    
    Args:
        pdb_file: Path to PDB file
        
    Returns:
        tuple: (antibody_ranges, antigen_ranges)
               where each ranges is a list of tuples (start, end)
    """
    
    # List of residues that are NOT protein (water, ions, solvent)
    solvent = {'WAT', 'HOH', 'TIP3', 'TIP4', 'SOL', 'Na+', 'Cl-', 'NA', 'CL', 'K+', 'MG2', 'Mg2', 'MG'}
    
    fragments = []  # List of (start, end, first_resname)
    current_fragment_start = None
    current_fragment_end = None
    current_fragment_first_resname = None
    prev_line_was_ter = False
    
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                resid = int(line[22:26].strip())
                resname = line[17:20].strip()
                
                # Ignore water and ions
                if resname in solvent:
                    continue
                
                # Start of new fragment
                if current_fragment_start is None or prev_line_was_ter:
                    # Save previous fragment if it exists
                    if current_fragment_start is not None:
                        fragments.append((
                            current_fragment_start,
                            current_fragment_end,
                            current_fragment_first_resname
                        ))
                    
                    # Start new fragment
                    current_fragment_start = resid
                    current_fragment_end = resid
                    current_fragment_first_resname = resname
                    prev_line_was_ter = False
                else:
                    # Continue current fragment
                    current_fragment_end = resid
            
            elif line.startswith('TER'):
                prev_line_was_ter = True
    
    # Save last fragment
    if current_fragment_start is not None:
        fragments.append((
            current_fragment_start,
            current_fragment_end,
            current_fragment_first_resname
        ))
    
    # Classify fragments into antibody and antigen
    antibody_ranges = []
    antigen_ranges = []
    
    for start, end, first_res in fragments:
        if first_res == 'ACE':
            antigen_ranges.append((start, end))
        else:
            antibody_ranges.append((start, end))
    
    return antibody_ranges, antigen_ranges

def format_ranges_for_cpptraj(ranges):
    """
    Converts list of ranges to cpptraj format.
    Example: [(1,100), (150,200)] -> "1-100,150-200"
    """
    if not ranges:
        return ""
    return ",".join([f"{start}-{end}" for start, end in ranges])

def generate_compaction_script(prmtop_in, traj_files, prmtop_out, traj_out, strip_mask=":WAT,Na+,Cl-"):
    """
    Generates compaction script to remove water and ions.
    
    Args:
        prmtop_in: Original topology (with water)
        traj_files: List of trajectory files or single pattern
        prmtop_out: Output topology (without water)
        traj_out: Output trajectory (without water)
        strip_mask: Mask for removal (default: water and ions)
    """
    # Generate multiple trajin lines if there are multiple files
    if isinstance(traj_files, list):
        trajin_lines = "\n".join([f"trajin {traj}" for traj in traj_files])
    else:
        trajin_lines = f"trajin {traj_files}"
    
    script = f"""# Step 1: Trajectory compaction
parm {prmtop_in}
{trajin_lines}
autoimage
strip {strip_mask} parmout {prmtop_out}
trajout {traj_out}
run
clear all
"""
    return script

def generate_rmsd_script(prmtop, traj, reference, ab_ranges, ag_ranges, output_rmsd, cutoff=8.0):
    """
    Generates RMSD analysis script.
    
    Args:
        prmtop: Compacted topology
        traj: Compacted trajectory
        reference: Reference file
        ab_ranges: List of tuples (start, end) for antibody
        ag_ranges: List of tuples (start, end) for antigen
        output_rmsd: Output file
        cutoff: Interface cutoff
    """
    ab_mask = format_ranges_for_cpptraj(ab_ranges)
    ag_mask = format_ranges_for_cpptraj(ag_ranges)
    
    script = f"""
# Step 2: iRMSD calculation
parm {prmtop}
trajin {traj}
autoimage
reference {reference}
rmsd :{ag_mask}<:{cutoff}&:{ab_mask}@CA|:{ab_mask}<:{cutoff}&:{ag_mask}@CA reference :{ag_mask}<:{cutoff}&:{ab_mask}@CA|:{ab_mask}<:{cutoff}&:{ag_mask}@CA out {output_rmsd}
run
"""
    return script

def main():
    parser = argparse.ArgumentParser(
        description='Compact trajectories and automatically calculate iRMSD',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Complete usage example:
  python script.py system.pdb \\
    --prmtop-in HMR_system.prmtop \\
    --traj-in rep/trajectory/mdProd310.nc rep/trajectory/mdProd330.nc \\
    --prmtop-out rep0_noWat_hMD.prmtop \\
    --traj-out rep0_noWat_hMD.nc \\
    --reference nowat_system.rst7 \\
    --output iRMSD_r0.dat \\
    --run

Or using wildcard in bash (use quotes):
  python script.py system.pdb \\
    --prmtop-in HMR_system.prmtop \\
    --traj-in "rep/trajectory/mdProd*.nc" \\
    --prmtop-out rep0_noWat_hMD.prmtop \\
    --traj-out rep0_noWat_hMD.nc \\
    --reference nowat_system.rst7 \\
    --output iRMSD_r0.dat \\
    --run

Example just to see residues:
  python script.py system.pdb --info

Example without executing (only generate scripts):
  python script.py system.pdb [options] --save cpptraj_full.in
        """
    )
    
    # Required arguments
    parser.add_argument('pdb', help='Input PDB file')
    
    # Arguments for compaction
    parser.add_argument('--prmtop-in', help='Original topology (with water)')
    parser.add_argument('--traj-in', nargs='+', help='Input trajectories (accepts multiple files or wildcard pattern)')
    parser.add_argument('--prmtop-out', help='Output topology (without water)')
    parser.add_argument('--traj-out', help='Output trajectory (without water)')
    parser.add_argument('--strip', default=':WAT,Na+,Cl-',
                        help='Mask for removal (default: :WAT,Na+,Cl-)')
    
    # Arguments for RMSD analysis
    parser.add_argument('-r', '--reference', help='Reference file (.rst7)')
    parser.add_argument('-o', '--output', help='RMSD output file name')
    parser.add_argument('-c', '--cutoff', type=float, default=8.0, 
                        help='Distance cutoff for interface (default: 8.0 Å)')
    
    # Control arguments
    parser.add_argument('--info', action='store_true', 
                        help='Only show residues without generating script')
    parser.add_argument('-s', '--save', help='Save complete cpptraj script to file')
    parser.add_argument('--run', action='store_true',
                        help='Execute cpptraj automatically')
    parser.add_argument('--compact-only', action='store_true',
                        help='Execute only compaction (no RMSD analysis)')
    parser.add_argument('--rmsd-only', action='store_true',
                        help='Execute only RMSD analysis (assumes trajectory already compacted)')
    
    args = parser.parse_args()
    
    # Extract and show residues
    antibody_ranges, antigen_ranges = extract_residue_ranges(args.pdb)
    
    print("=" * 60)
    print(f"Fragments detected in {args.pdb}")
    print("=" * 60)
    print(f"\n🔵 ANTIBODY ({len(antibody_ranges)} fragment(s)):")
    print("   (fragments that DON'T start with ACE)")
    for i, (start, end) in enumerate(antibody_ranges, 1):
        print(f"   Fragment {i}: residues {start}-{end}")
    
    print(f"\n🔴 ANTIGEN ({len(antigen_ranges)} fragment(s)):")
    print("   (fragments that start with ACE)")
    for i, (start, end) in enumerate(antigen_ranges, 1):
        print(f"   Fragment {i}: residues {start}-{end}")
    
    # Format ranges for cpptraj
    ab_mask = format_ranges_for_cpptraj(antibody_ranges)
    ag_mask = format_ranges_for_cpptraj(antigen_ranges)
    
    print(f"\n📋 Masks for cpptraj:")
    print(f"   Antibody: :{ab_mask}")
    print(f"   Antigen:  :{ag_mask}")
    print("=" * 60)
    print()
    
    # Validations
    if len(antibody_ranges) != 2:
        print(f"⚠️  WARNING: Expected 2 antibody fragments, found {len(antibody_ranges)}")
        print("   Check if capping was done correctly")
        print()
    
    if len(antigen_ranges) == 0:
        print("❌ ERROR: No antigen fragments detected!")
        print("   Check if antigen has ACE at N-terminals")
        return
    
    if args.info:
        return
    
    # Check necessary parameters
    if not args.rmsd_only:
        if not all([args.prmtop_in, args.traj_in, args.prmtop_out, args.traj_out]):
            print("❌ For compaction, provide:")
            print("   --prmtop-in, --traj-in, --prmtop-out, --traj-out")
            print()
            print("Use --info to only see residues")
            print("Use --rmsd-only if trajectory is already compacted")
            return
    
    if not args.compact_only:
        if not all([args.reference, args.output]):
            print("❌ For RMSD analysis, provide:")
            print("   --reference, --output")
            print()
            print("Use --compact-only to only compact")
            return
        
        if args.rmsd_only and not args.prmtop_out:
            print("❌ For --rmsd-only, provide --prmtop-out (compacted topology)")
            return
    
    # Process input trajectories
    if args.traj_in:
        # If it's a single element list with wildcard, expand via glob
        if len(args.traj_in) == 1 and ('*' in args.traj_in[0] or '?' in args.traj_in[0]):
            import glob
            expanded = sorted(glob.glob(args.traj_in[0]))
            if expanded:
                args.traj_in = expanded
                print(f"🔍 Found {len(expanded)} trajectories:")
                for traj in expanded:
                    print(f"   - {traj}")
                print()
            else:
                print(f"⚠️  No files found for pattern: {args.traj_in[0]}")
                return
    
    # Generate complete script
    full_script = ""
    
    if not args.rmsd_only:
        print("📦 Generating compaction script...")
        compact_script = generate_compaction_script(
            args.prmtop_in, args.traj_in, args.prmtop_out, 
            args.traj_out, args.strip
        )
        full_script += compact_script
        print("✓ Compaction script generated")
        print()
    
    if not args.compact_only:
        print("📊 Generating RMSD analysis script...")
        # Define which topology and trajectory to use
        prmtop_rmsd = args.prmtop_out if not args.rmsd_only else args.prmtop_out
        traj_rmsd = args.traj_out if not args.rmsd_only else args.traj_in
        
        rmsd_script = generate_rmsd_script(
            prmtop_rmsd, traj_rmsd, args.reference,
            antibody_ranges, antigen_ranges,
            args.output, args.cutoff
        )
        full_script += rmsd_script
        print("✓ RMSD analysis script generated")
        print()
    
    # Show generated script
    print("=" * 60)
    print("Complete cpptraj script")
    print("=" * 60)
    print(full_script)
    print("=" * 60)
    print()
    
    # Save to file
    script_file = args.save if args.save else 'cpptraj_full.in'
    with open(script_file, 'w') as f:
        f.write(full_script)
    print(f"💾 Script saved to: {script_file}")
    print()
    
    # Execute cpptraj if requested
    if args.run:
        print("=" * 60)
        print("Executing cpptraj")
        print("=" * 60)
        try:
            result = subprocess.run(['cpptraj', '-i', script_file], 
                                  capture_output=True, text=True)
            print(result.stdout)
            if result.stderr:
                print("STDERR:", result.stderr)
            
            # Check if output files were created
            files_created = []
            files_missing = []
            
            if not args.compact_only and args.output:
                if os.path.exists(args.output):
                    files_created.append(args.output)
                else:
                    files_missing.append(args.output)
            
            if not args.rmsd_only:
                if os.path.exists(args.prmtop_out):
                    files_created.append(args.prmtop_out)
                else:
                    files_missing.append(args.prmtop_out)
                    
                if os.path.exists(args.traj_out):
                    files_created.append(args.traj_out)
                else:
                    files_missing.append(args.traj_out)
            
            if result.returncode == 0 and not files_missing:
                print()
                print("=" * 60)
                print("✅ Analysis completed successfully!")
                print("=" * 60)
                if files_created:
                    print("📁 Files created:")
                    for f in files_created:
                        size = os.path.getsize(f) / (1024*1024)  # MB
                        print(f"   ✓ {f} ({size:.2f} MB)")
            else:
                print()
                print("=" * 60)
                if files_missing:
                    print("⚠️  Execution completed but some files were not created:")
                    print("=" * 60)
                    for f in files_missing:
                        print(f"   ✗ {f}")
                    print()
                    print("💡 Possible causes:")
                    print("   - Error in cpptraj (check output above)")
                    print("   - Incorrect file path")
                    print("   - Invalid residue selection")
                    print("   - Problem with reference file")
                else:
                    print(f"❌ Execution error (code: {result.returncode})")
                print("=" * 60)
        except FileNotFoundError:
            print("❌ Error: cpptraj not found in PATH")
            print(f"Execute manually: cpptraj -i {script_file}")
    else:
        print(f"To execute: cpptraj -i {script_file}")

if __name__ == "__main__":
    main()
