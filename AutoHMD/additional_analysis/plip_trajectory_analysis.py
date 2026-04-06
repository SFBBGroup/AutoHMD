#!/usr/bin/env python3
"""
PLIP Trajectory Analysis - Unified Script
Analyzes protein-protein interactions across MD trajectory frames using PLIP.

MODES:
  1. Config file (auto-detect):   python plip_trajectory_analysis.py config.yaml
  2. Config file (manual chains): python plip_trajectory_analysis.py config.yaml
                                  (set auto_detect_chains: false + residue_ranges in config)
  3. Command line (auto-detect):  python plip_trajectory_analysis.py \\
                                    --topology topology.prmtop \\
                                    --trajectories traj1.nc traj2.nc \\
                                    --output results/ --interval 32

Dependencies:
    - cpptraj (AmberTools)
    - pdb-tools
    - PLIP
    - Python packages: pyyaml, pandas, matplotlib, seaborn, numpy, tqdm
"""

import os
import re
import sys
import yaml
import shutil
import subprocess
import tempfile
import argparse
from pathlib import Path
import xml.etree.ElementTree as ET
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict
from tqdm import tqdm


# ============================================================
# AUTO-DETECT CHAINS
# ============================================================

def extract_reference_pdb_from_topology(topology_file, output_pdb, trajectory_file=None):
    """Extract a single-frame PDB from topology (and optionally first frame of trajectory)."""
    print(f"\nExtracting reference PDB from topology...")

    if trajectory_file is None:
        script = f"parm {topology_file}\ntrajout {output_pdb} pdb\nrun\nquit\n"
    else:
        script = (
            f"parm {topology_file}\n"
            f"trajin {trajectory_file} 1 1\n"
            "strip :WAT,Na+,Cl-\n"
            f"trajout {output_pdb} pdb\n"
            "run\nquit\n"
        )

    script_file = Path(output_pdb).parent / "extract_ref.cpptraj"
    script_file.write_text(script)

    result = subprocess.run(["cpptraj", "-i", str(script_file)], capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"Failed to extract PDB from topology:\n{result.stderr}")
    if not Path(output_pdb).exists():
        raise RuntimeError(f"Reference PDB not created: {output_pdb}")

    print(f"Reference PDB extracted: {output_pdb}")
    return output_pdb


def detect_chains_from_pdb(pdb_file):
    """
    Auto-detect chain assignments from a PDB file.

    Logic:
      - Fragments are delimited by NME residues.
      - A fragment starting with ACE  → antigen (chain B)
      - Fragments NOT starting with ACE → antibody (first = H, second = L)
    """
    print(f"\n{'='*60}")
    print("AUTOMATIC CHAIN DETECTION")
    print(f"{'='*60}")
    print(f"Reading PDB file: {pdb_file}")

    residues = []
    current_res = None
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith(('ATOM', 'HETATM')):
                res_name = line[17:20].strip()
                res_num  = int(line[22:26].strip())
                if current_res != res_num:
                    residues.append({'num': res_num, 'name': res_name})
                    current_res = res_num

    if not residues:
        raise RuntimeError("No residues found in PDB file")

    print(f"Total residues found: {len(residues)}")
    print(f"First residue: {residues[0]['name']}{residues[0]['num']}")
    print(f"Last residue:  {residues[-1]['name']}{residues[-1]['num']}")

    # Split into fragments at NME residues
    fragments = []
    current_fragment = []
    for res in residues:
        current_fragment.append(res)
        if res['name'] == 'NME':
            fragments.append(current_fragment)
            current_fragment = []
    if current_fragment:
        fragments.append(current_fragment)

    print(f"\nDetected {len(fragments)} fragments:")

    antibody_fragments = []
    antigen_fragments  = []

    for i, frag in enumerate(fragments):
        starts_with_ace = frag[0]['name'] == 'ACE'
        ends_with_nme   = frag[-1]['name'] == 'NME'
        start_res = frag[0]['num']
        end_res   = frag[-1]['num']

        label = (f"  Fragment {i+1}: residues {start_res}-{end_res} "
                 f"({len(frag)} residues, "
                 f"{'ACE->' if starts_with_ace else frag[0]['name']}"
                 f"{'...->NME' if ends_with_nme else ''})")
        print(label)

        if starts_with_ace:
            antigen_fragments.append((start_res, end_res))
        else:
            antibody_fragments.append((start_res, end_res))

    if len(antibody_fragments) < 2:
        raise RuntimeError(
            f"Expected at least 2 antibody fragments (H and L), found {len(antibody_fragments)}"
        )
    if not antigen_fragments:
        raise RuntimeError("No antigen fragments detected (no ACE residues found)")

    heavy_chain = antibody_fragments[0]
    light_chain = antibody_fragments[1]
    antigen_range = (
        min(s for s, _ in antigen_fragments),
        max(e for _, e in antigen_fragments),
    )

    print(f"\n{'='*60}")
    print("CHAIN ASSIGNMENT:")
    print(f"{'='*60}")
    print(f"  Heavy Chain (H): residues {heavy_chain[0]}-{heavy_chain[1]}")
    print(f"  Light Chain (L): residues {light_chain[0]}-{light_chain[1]}")
    print(f"  Antigen     (B): residues {antigen_range[0]}-{antigen_range[1]}")

    if len(antibody_fragments) > 2:
        print(f"\n  WARNING: Found {len(antibody_fragments)} antibody fragments; using first 2")

    return {
        'H': list(heavy_chain),
        'L': list(light_chain),
        'B': list(antigen_range),
    }


# ============================================================
# CONFIGURATION
# ============================================================

class Config:
    """
    Unified configuration class.

    Accepts:
      - config_file (str/Path) : YAML config with either manual or auto-detect settings
      - args (Namespace)       : parsed CLI arguments (always uses auto-detect)
    """

    def __init__(self, config_file=None, args=None):
        if args is not None:
            self._from_args(args)
        elif config_file is not None:
            self._from_file(config_file)
        else:
            raise ValueError("Provide either config_file or args")

    # ----------------------------------------------------------
    def _from_file(self, config_file):
        with open(config_file, 'r') as f:
            data = yaml.safe_load(f)

        traj_block = data['trajectory']
        self.topology = traj_block['topology']

        if 'trajectories' in traj_block:
            self.trajectories = traj_block['trajectories']
        elif 'trajectory' in traj_block:
            self.trajectories = [traj_block['trajectory']]
        else:
            raise ValueError("Config must specify 'trajectory' or 'trajectories'")

        self.reference = traj_block.get('reference', None)

        frames = data.get('frames', {})
        self.start_frame = frames.get('start', 1)
        self.end_frame   = frames.get('end', -1)
        self.interval    = frames.get('interval', 10)

        # --- Chain assignment mode ---
        # Priority: auto_detect_chains flag > residue_ranges > reference_for_chains
        self.auto_detect_chains  = data.get('auto_detect_chains', False)
        self.reference_for_chains = data.get('reference_for_chains', None)
        self.residue_ranges       = data.get('residue_ranges', {})

        # PLIP
        self.target_chain = data.get('plip', {}).get('target_chain', 'B')

        # Output
        self.output_dir    = Path(data['output']['directory'])
        self.output_prefix = data['output']['prefix']

        # PDB processing
        self.chains_to_keep   = data.get('chains_to_keep', ['H', 'L', 'B'])
        self.reres_starts     = data.get('reres_starts', {'H': 1, 'L': 1, 'B': 1})
        self.fixinsert_chains = set(data.get('fixinsert_chains', ['H', 'L']))

        self.replica_size = data.get('replica_size', 175)

    # ----------------------------------------------------------
    def _from_args(self, args):
        self.topology     = args.topology
        self.trajectories = args.trajectories
        self.reference    = None

        self.start_frame = args.start_frame
        self.end_frame   = args.end_frame
        self.interval    = args.interval

        # CLI mode always uses auto-detect
        self.auto_detect_chains   = True
        self.reference_for_chains = None
        self.residue_ranges       = {}

        self.target_chain  = args.target_chain
        self.output_dir    = Path(args.output)
        self.output_prefix = args.prefix

        self.chains_to_keep   = ['H', 'L', 'B']
        self.reres_starts     = {'H': 1, 'L': 1, 'B': 1}
        self.fixinsert_chains = {'H', 'L'}

        self.replica_size = args.replica_size

    # ----------------------------------------------------------
    def resolve_chain_mode(self, work_dir):
        """
        Decide which chain-assignment strategy to use and populate
        self.residue_ranges if auto-detect is requested.
        """
        print(f"\n{'='*60}")
        print("CHAIN ASSIGNMENT MODE")
        print(f"{'='*60}")

        if self.auto_detect_chains and not self.residue_ranges:
            print("  Mode: AUTO-DETECT (from topology/trajectory)")
            ref_pdb = work_dir / "reference_for_detection.pdb"
            extract_reference_pdb_from_topology(
                self.topology, ref_pdb, self.trajectories[0]
            )
            self.residue_ranges = detect_chains_from_pdb(ref_pdb)
            self.reference_for_chains = None
            print("\n  Auto-detected residue ranges:")
            for chain, rng in self.residue_ranges.items():
                print(f"    Chain {chain}: {rng[0]}-{rng[1]}")

        elif self.residue_ranges:
            print("  Mode: MANUAL (residue ranges from config)")
            for chain, rng in self.residue_ranges.items():
                print(f"    Chain {chain}: {rng}")

        elif self.reference_for_chains:
            print(f"  Mode: REFERENCE PDB  ({self.reference_for_chains})")

        else:
            print("  WARNING: No chain assignment method configured!")
            print("  PDB will use whatever chain IDs cpptraj generated.")


# ============================================================
# CPPTRAJ INTERFACE
# ============================================================

def generate_cpptraj_script(config, frame_dir):
    end_str = str(config.end_frame) if config.end_frame > 0 else "last"
    script = f"# PLIP trajectory analysis - frame extraction\nparm {config.topology}\n"
    for traj in config.trajectories:
        script += f"trajin {traj} {config.start_frame} {end_str} {config.interval}\n"
    if config.reference:
        script += f"reference {config.reference}\nrms reference\n"
    script += f"autoimage\nstrip :WAT,Na+,Cl-\ntrajout {frame_dir}/frame.pdb pdb multi\nrun\nquit\n"
    return script


def extract_frames(config, work_dir):
    frame_dir = work_dir / "frames"
    frame_dir.mkdir(exist_ok=True)

    script_file = work_dir / "extract_frames.cpptraj"
    script_file.write_text(generate_cpptraj_script(config, frame_dir))

    print(f"\n{'='*60}")
    print("EXTRACTING FRAMES WITH CPPTRAJ")
    print(f"{'='*60}")
    print(f"Script: {script_file}")

    result = subprocess.run(["cpptraj", "-i", str(script_file)], capture_output=True, text=True)
    if result.returncode != 0:
        print("CPPTRAJ STDERR:")
        print(result.stderr)
        raise RuntimeError(f"cpptraj failed with return code {result.returncode}")

    frames = []
    for pattern in ["frame.pdb.*", "frame_.pdb.*", "frame_*.pdb", "*.pdb"]:
        frames = [f for f in frame_dir.glob(pattern) if f.name != "frame.pdb"]
        if frames:
            break

    if not frames:
        all_files = list(frame_dir.glob("*"))
        print(f"\nFiles in {frame_dir}:")
        for f in all_files:
            print(f"  {f.name}")
        raise RuntimeError(f"No PDB frames found in {frame_dir}")

    def _frame_num(path):
        m = re.search(r'(\d+)', path.name)
        return int(m.group(1)) if m else 0

    frames = sorted(frames, key=_frame_num)
    print(f"Extracted {len(frames)} frames")
    return frames


# ============================================================
# CHAIN ASSIGNMENT
# ============================================================

def assign_chains_from_reference(pdb_file, reference_pdb):
    with open(reference_pdb, 'r') as f:
        ref_lines = [l for l in f if l.startswith(('ATOM', 'HETATM'))]
    output_lines = []
    with open(pdb_file, 'r') as f:
        atom_idx = 0
        for line in f:
            if line.startswith(('ATOM', 'HETATM')) and atom_idx < len(ref_lines):
                output_lines.append(line[:21] + ref_lines[atom_idx][21] + line[22:])
                atom_idx += 1
            else:
                output_lines.append(line)
    with open(pdb_file, 'w') as f:
        f.writelines(output_lines)


def assign_chains_from_ranges(pdb_file, residue_ranges):
    output_lines = []
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith(('ATOM', 'HETATM')):
                try:
                    resnum = int(line[22:26].strip())
                    chain  = ' '
                    for chain_id, rng in residue_ranges.items():
                        start, end = rng[0], rng[1]
                        if start <= resnum <= end:
                            chain = chain_id
                            break
                    output_lines.append(line[:21] + chain + line[22:])
                except (ValueError, IndexError):
                    output_lines.append(line)
            else:
                output_lines.append(line)
    with open(pdb_file, 'w') as f:
        f.writelines(output_lines)


def fix_chain_assignments(pdb_file, config):
    if config.reference_for_chains and Path(config.reference_for_chains).exists():
        assign_chains_from_reference(pdb_file, config.reference_for_chains)
    elif config.residue_ranges:
        assign_chains_from_ranges(pdb_file, config.residue_ranges)
    else:
        print("  WARNING: No chain assignment method configured; using cpptraj chain IDs")


# ============================================================
# PDB PROCESSING (pdb-tools pipeline)
# ============================================================

def check_pdb_tools():
    required = ["pdb_keepcoord", "pdb_selchain", "pdb_reres", "pdb_fixinsert"]
    missing = [t for t in required if shutil.which(t) is None]
    if missing:
        raise RuntimeError(f"Missing required pdb-tools: {', '.join(missing)}")


def _run_cmd(cmd, input_bytes=None):
    proc = subprocess.run(cmd, input=input_bytes, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if proc.returncode != 0:
        raise RuntimeError(
            f"Command failed: {' '.join(cmd)}\n{proc.stderr.decode(errors='replace')}"
        )
    return proc.stdout


def _process_chain(in_pdb, chain_id, tmpdir, reres_start, fixinsert):
    tmp_keep  = tmpdir / f"{chain_id}_keep.pdb"
    tmp_sel   = tmpdir / f"{chain_id}_sel.pdb"
    tmp_reres = tmpdir / f"{chain_id}_reres.pdb"
    tmp_final = tmpdir / f"{chain_id}_final.pdb"

    stdout  = _run_cmd(["pdb_keepcoord", str(in_pdb)])
    cleaned = b"\n".join(l for l in stdout.splitlines()
                         if l.startswith((b"ATOM", b"HETATM"))) + b"\n"
    tmp_keep.write_bytes(cleaned)

    tmp_sel.write_bytes(_run_cmd(["pdb_selchain", f"-{chain_id}", str(tmp_keep)]))
    tmp_reres.write_bytes(_run_cmd(["pdb_reres", f"-{reres_start}", str(tmp_sel)]))

    if fixinsert:
        tmp_final.write_bytes(_run_cmd(["pdb_fixinsert", str(tmp_reres)]))
    else:
        tmp_final.write_bytes(tmp_reres.read_bytes())

    return tmp_final


def process_pdb_file(input_pdb, output_pdb, config):
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        processed = []
        for ch in config.chains_to_keep:
            try:
                p = _process_chain(
                    input_pdb, ch, td,
                    config.reres_starts.get(ch, 1),
                    ch in config.fixinsert_chains
                )
                processed.append((ch, p))
            except RuntimeError:
                continue
        if not processed:
            raise RuntimeError(f"No chains processed for {input_pdb}")

        lines = []
        for ch, path in processed:
            txt = path.read_text().splitlines()
            while txt and txt[-1].strip().upper().startswith(("END", "ENDMDL")):
                txt.pop()
            if txt and not txt[-1].startswith("TER"):
                txt.append("TER")
            lines.extend(txt)
        lines.append("END")
        Path(output_pdb).write_bytes(("\n".join(lines) + "\n").encode())


# ============================================================
# PLIP INTERFACE
# ============================================================

def run_plip(pdb_file, target_chain, output_dir):
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True, parents=True)
    result = subprocess.run(
        ["plip", "-f", str(pdb_file), "--inter", target_chain, "-vx", "-o", str(output_dir)],
        capture_output=True, text=True
    )
    if result.returncode != 0:
        raise RuntimeError(f"PLIP failed on {pdb_file}:\n{result.stderr}")
    xml_files = list(output_dir.glob("**/*report.xml"))
    if not xml_files:
        raise RuntimeError(f"No report.xml found in {output_dir}")
    return xml_files[0]


# ============================================================
# XML PARSING
# ============================================================

def _safe_float(val, default=np.nan):
    try:
        return float(val) if val is not None else default
    except (TypeError, ValueError):
        return default


def _find_text(elem, tag, default=None):
    return elem.findtext(tag, default) if elem is not None else default


def parse_plip_xml(xml_file):
    tree = ET.parse(xml_file)
    root = tree.getroot()

    interactions = {k: [] for k in
                    ['hydrophobic', 'hydrogen_bonds', 'salt_bridges',
                     'pi_stacks', 'pi_cation', 'halogen_bonds', 'metal_complexes']}

    section = root.find('.//interactions')
    if section is None:
        return interactions

    parsers = {
        'hydrophobic':    ('hydrophobic_interactions',  'hydrophobic_interaction',  'Hydrophobic'),
        'hydrogen_bonds': ('hydrogen_bonds',             'hydrogen_bond',            'Hydrogen Bond'),
        'salt_bridges':   ('salt_bridges',               'salt_bridge',              'Salt Bridge'),
        'pi_stacks':      ('pi_stacks',                  'pi_stack',                 'Pi-Pi Stacking'),
        'pi_cation':      ('pi_cation_interactions',     'pi_cation_interaction',    'Pi-Cation'),
    }

    for key, (parent_tag, child_tag, itype) in parsers.items():
        parent = section.find(parent_tag)
        if parent is None:
            continue
        for item in parent.findall(child_tag):
            try:
                interactions[key].append({
                    'prot_res':   f"{_find_text(item,'restype','').strip()}{_find_text(item,'resnr','').strip()}",
                    'prot_chain': _find_text(item, 'reschain'),
                    'lig_res':    f"{_find_text(item,'restype_lig','').strip()}{_find_text(item,'resnr_lig','').strip()}",
                    'lig_chain':  _find_text(item, 'reschain_lig'),
                    'distance':   _safe_float(_find_text(item, 'dist')),
                    'type':       itype,
                })
            except Exception:
                continue

    return interactions


# ============================================================
# AGGREGATION & PREVALENCE
# ============================================================

def aggregate_interactions(all_frame_data):
    all_interactions = []
    for frame_num, interactions in all_frame_data.items():
        for int_list in interactions.values():
            for interaction in int_list:
                interaction['frame'] = frame_num
                all_interactions.append(interaction)
    if not all_interactions:
        return None
    df = pd.DataFrame(all_interactions)
    df['interaction_id'] = df.apply(
        lambda r: f"{r['prot_res']}_{r['prot_chain']}--{r['lig_res']}_{r['lig_chain']}_{r['type']}",
        axis=1
    )
    return df


def calculate_prevalence(df, total_frames):
    if df is None or df.empty:
        return None
    prev = df.groupby('interaction_id').agg(
        frame_count=('frame', 'nunique'),
        type=('type', 'first'),
        prot_res=('prot_res', 'first'),
        prot_chain=('prot_chain', 'first'),
        lig_res=('lig_res', 'first'),
        lig_chain=('lig_chain', 'first'),
        avg_distance=('distance', 'mean'),
    ).reset_index()
    prev['prevalence_pct'] = (prev['frame_count'] / total_frames) * 100
    total_occ = df.groupby('interaction_id').size()
    prev['total_occurrences'] = prev['interaction_id'].map(total_occ)
    prev['avg_per_frame'] = prev['total_occurrences'] / prev['frame_count']
    return prev.sort_values('prevalence_pct', ascending=False)


# ============================================================
# VISUALIZATION
# ============================================================

def plot_prevalence_overview(prevalence_df, output_file):
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))

    ax = axes[0, 0]
    top_20 = prevalence_df.head(20)
    ax.barh(range(len(top_20)), top_20['prevalence_pct'],
            color=sns.color_palette("husl", len(top_20)))
    ax.set_yticks(range(len(top_20)))
    ax.set_yticklabels([f"{r['prot_res']}--{r['lig_res']}" for _, r in top_20.iterrows()],
                       fontsize=8)
    ax.set_xlabel('Prevalence (%)')
    ax.set_title('Top 20 Interactions by Prevalence', fontweight='bold')
    ax.invert_yaxis()

    ax = axes[0, 1]
    type_counts = prevalence_df['type'].value_counts()
    ax.pie(type_counts.values, labels=type_counts.index, autopct='%1.1f%%',
           colors=sns.color_palette("Set2", len(type_counts)), startangle=90)
    ax.set_title('Interaction Type Distribution', fontweight='bold')

    ax = axes[1, 0]
    ax.hist(prevalence_df['prevalence_pct'], bins=20, color='skyblue', edgecolor='black')
    ax.set_xlabel('Prevalence (%)')
    ax.set_ylabel('Number of Interactions')
    ax.set_title('Distribution of Interaction Prevalence', fontweight='bold')
    median = prevalence_df['prevalence_pct'].median()
    ax.axvline(median, color='red', linestyle='--', label=f'Median: {median:.1f}%')
    ax.legend()

    ax = axes[1, 1]
    for itype in prevalence_df['type'].unique():
        sub = prevalence_df[prevalence_df['type'] == itype]
        ax.scatter(sub['prevalence_pct'], sub['avg_distance'], label=itype, alpha=0.6, s=50)
    ax.set_xlabel('Prevalence (%)')
    ax.set_ylabel('Average Distance (Å)')
    ax.set_title('Distance vs Prevalence by Interaction Type', fontweight='bold')
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_file}")


def plot_interaction_timeline(df, prevalence_df, output_file, top_n=15):
    top_ids = prevalence_df.head(top_n)['interaction_id'].tolist()
    df_top  = df[df['interaction_id'].isin(top_ids)]
    frames  = sorted(df['frame'].unique())

    presence_matrix = np.zeros((len(top_ids), len(frames)))
    for i, int_id in enumerate(top_ids):
        for j, frame in enumerate(frames):
            if int_id in df_top[df_top['frame'] == frame]['interaction_id'].values:
                presence_matrix[i, j] = 1

    fig, ax = plt.subplots(figsize=(16, 10))
    im = ax.imshow(presence_matrix, aspect='auto', cmap='YlOrRd', interpolation='nearest')
    ax.set_yticks(range(len(top_ids)))
    y_labels = []
    for int_id in top_ids:
        row = prevalence_df[prevalence_df['interaction_id'] == int_id].iloc[0]
        y_labels.append(f"{row['prot_res']}--{row['lig_res']} ({row['type'][:4]})")
    ax.set_yticklabels(y_labels, fontsize=9)
    ax.set_xlabel('Frame Number', fontweight='bold')
    ax.set_ylabel('Interaction', fontweight='bold')
    ax.set_title(f'Top {top_n} Interactions Timeline', fontweight='bold', pad=20)
    plt.colorbar(im, ax=ax).set_label('Present', rotation=270, labelpad=15)
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_file}")


def plot_residue_heatmap(prevalence_df, output_file):
    residue_pairs = defaultdict(int)
    for _, row in prevalence_df.iterrows():
        residue_pairs[(row['prot_res'], row['lig_res'])] += row['prevalence_pct']

    if not residue_pairs:
        print("No residue pairs to plot")
        return

    prot_res = sorted(set(k[0] for k in residue_pairs))
    lig_res  = sorted(set(k[1] for k in residue_pairs))
    matrix = np.zeros((len(prot_res), len(lig_res)))
    for i, pr in enumerate(prot_res):
        for j, lr in enumerate(lig_res):
            matrix[i, j] = residue_pairs.get((pr, lr), 0)

    fig, ax = plt.subplots(figsize=(max(12, len(lig_res)*0.4), max(10, len(prot_res)*0.3)))
    sns.heatmap(matrix, xticklabels=lig_res, yticklabels=prot_res,
                cmap='YlOrRd', annot=False, cbar_kws={'label': 'Total Prevalence (%)'})
    ax.set_xlabel('Antigen Residues', fontweight='bold')
    ax.set_ylabel('Antibody Residues', fontweight='bold')
    ax.set_title('Residue-Residue Interaction Prevalence', fontweight='bold', pad=20)
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_file}")


def plot_interaction_map(prevalence_df, output_file, min_prevalence=10):
    df_f = prevalence_df[prevalence_df['prevalence_pct'] >= min_prevalence].copy()
    if df_f.empty:
        print(f"No interactions with prevalence >= {min_prevalence}%")
        return

    int_type_colors = {
        'Hydrogen Bond': '#1f77b4', 'Hydrophobic': '#ff7f0e',
        'Salt Bridge': '#d62728',   'Pi-Pi Stacking': '#9467bd',
        'Pi-Cation': '#2ca02c',     'Halogen Bond': '#8c564b',
        'Metal Complex': '#e377c2'
    }

    fig = plt.figure(figsize=(20, 10))
    ax1 = plt.subplot(1, 2, 1)

    ab_res = sorted(df_f['prot_res'].unique())
    ag_res = sorted(df_f['lig_res'].unique())
    ab_pos = {r: i for i, r in enumerate(ab_res)}
    ag_pos = {r: i for i, r in enumerate(ag_res)}

    for _, row in df_f.iterrows():
        color = int_type_colors.get(row['type'], 'gray')
        lw    = 0.5 + (row['prevalence_pct'] / 100) * 4
        alpha = 0.3 + (row['prevalence_pct'] / 100) * 0.6
        ax1.plot([0, 1], [ab_pos[row['prot_res']], ag_pos[row['lig_res']]],
                 color=color, linewidth=lw, alpha=alpha, zorder=1)

    max_y = max(len(ab_res), len(ag_res))
    for res, pos in ab_pos.items():
        ax1.scatter(0, pos, s=100, c='skyblue', edgecolor='black', zorder=2)
        ax1.text(-0.05, pos, res, ha='right', va='center', fontsize=9, fontweight='bold')
    for res, pos in ag_pos.items():
        ax1.scatter(1, pos, s=100, c='lightcoral', edgecolor='black', zorder=2)
        ax1.text(1.05, pos, res, ha='left', va='center', fontsize=9, fontweight='bold')

    ax1.set_xlim(-0.3, 1.3)
    ax1.set_ylim(-1, max_y)
    ax1.axis('off')
    ax1.set_title(f'Interaction Map (prevalence ≥ {min_prevalence}%)', fontsize=14, fontweight='bold', pad=20)
    ax1.text(0, max_y + 0.5, 'Antibody', ha='center', fontsize=12, fontweight='bold')
    ax1.text(1, max_y + 0.5, 'Antigen',  ha='center', fontsize=12, fontweight='bold')
    legend_el = [plt.Line2D([0],[0], color=c, linewidth=3, label=t)
                 for t, c in int_type_colors.items() if t in df_f['type'].values]
    ax1.legend(handles=legend_el, loc='lower left', fontsize=10)

    ax2 = plt.subplot(1, 2, 2)
    ax2.axis('tight'); ax2.axis('off')
    table_data = []
    for i, (_, row) in enumerate(df_f.head(20).iterrows(), 1):
        table_data.append([f"{i}", row['prot_res'], row['lig_res'],
                           row['type'][:10], f"{row['prevalence_pct']:.1f}%",
                           f"{row['avg_distance']:.2f}Å"])
    t = ax2.table(cellText=table_data,
                  colLabels=['#', 'Ab Res', 'Ag Res', 'Type', 'Prev.', 'Dist.'],
                  cellLoc='center', loc='center',
                  colWidths=[0.08, 0.15, 0.15, 0.25, 0.15, 0.15])
    t.auto_set_font_size(False); t.set_fontsize(9); t.scale(1, 2)
    for j in range(6):
        t[(0, j)].set_facecolor('#4CAF50')
        t[(0, j)].set_text_props(weight='bold', color='white')
    for i in range(1, len(table_data)+1):
        for j in range(6):
            if i % 2 == 0:
                t[(i, j)].set_facecolor('#f0f0f0')
    ax2.set_title('Top 20 Interactions', fontsize=14, fontweight='bold', pad=20)

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_file}")


def plot_chain_interaction_matrix(prevalence_df, output_file):
    chain_ag = defaultdict(lambda: defaultdict(float))
    for _, row in prevalence_df.iterrows():
        chain_ag[row['prot_chain']][row['lig_res']] += row['prevalence_pct']

    if not chain_ag:
        print("No chain interaction data"); return

    ab_chains  = sorted(chain_ag.keys())
    ag_residues = sorted(set().union(*[chain_ag[c].keys() for c in ab_chains]))
    matrix = np.zeros((len(ab_chains), len(ag_residues)))
    for i, ch in enumerate(ab_chains):
        for j, res in enumerate(ag_residues):
            matrix[i, j] = chain_ag[ch][res]

    fig, ax = plt.subplots(figsize=(max(14, len(ag_residues)*0.4), 6))
    im = ax.imshow(matrix, aspect='auto', cmap='YlOrRd', interpolation='nearest')
    ax.set_xticks(range(len(ag_residues))); ax.set_xticklabels(ag_residues, rotation=45, ha='right')
    ax.set_yticks(range(len(ab_chains))); ax.set_yticklabels([f'Chain {c}' for c in ab_chains])
    ax.set_xlabel('Antigen Residues', fontweight='bold', fontsize=12)
    ax.set_ylabel('Antibody Chains',  fontweight='bold', fontsize=12)
    ax.set_title('Chain-Level Interaction Heatmap', fontweight='bold', fontsize=14, pad=20)
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('Total Prevalence (%)', rotation=270, labelpad=20, fontweight='bold')
    for i in range(len(ab_chains)):
        for j in range(len(ag_residues)):
            if matrix[i, j] > 0:
                tc = 'white' if matrix[i, j] > matrix.max()/2 else 'black'
                ax.text(j, i, f'{matrix[i, j]:.0f}', ha='center', va='center', color=tc, fontsize=8)
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_file}")


def plot_interaction_fingerprint(df, prevalence_df, output_file, top_n=15):
    top = prevalence_df.head(top_n)
    frames = sorted(df['frame'].unique())
    type_colors = {
        'Hydrogen Bond': '#1f77b4', 'Hydrophobic': '#ff7f0e',
        'Salt Bridge': '#d62728',   'Pi-Pi Stacking': '#9467bd',
        'Pi-Cation': '#2ca02c',
    }
    fig, axes = plt.subplots(top_n, 1, figsize=(16, top_n * 0.8), sharex=True)
    if top_n == 1:
        axes = [axes]

    for idx, (_, row) in enumerate(top.iterrows()):
        ax = axes[idx]
        int_data = df[df['interaction_id'] == row['interaction_id']]
        for frame in frames:
            fd = int_data[int_data['frame'] == frame]
            if not fd.empty:
                color = type_colors.get(fd.iloc[0]['type'], 'gray')
                ax.barh(0, 1, left=frame, height=0.8, color=color, edgecolor='none', alpha=0.8)
        ax.set_yticks([0])
        ax.set_yticklabels([f"{row['prot_res']}--{row['lig_res']}"], fontsize=9)
        ax.set_ylim(-0.5, 0.5)
        for spine in ['top', 'right', 'left']:
            ax.spines[spine].set_visible(False)
        ax.text(frames[-1] * 1.02, 0, f"{row['prevalence_pct']:.0f}%",
                va='center', fontsize=8, fontweight='bold')

    axes[-1].set_xlabel('Frame Number', fontweight='bold', fontsize=12)
    axes[0].set_title(f'Interaction Fingerprint - Top {top_n} Interactions',
                      fontsize=14, fontweight='bold', pad=20)
    legend_el = [plt.Rectangle((0,0),1,1, fc=c, label=t) for t, c in type_colors.items()]
    axes[0].legend(handles=legend_el, loc='upper right', ncol=len(type_colors),
                   fontsize=8, framealpha=0.9)
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_file}")


def plot_unique_interactions_over_time(df, output_file, replica_size=175):
    """
    Plot total unique interactions per frame and stacked interaction types over time.
    Vertical dashed lines mark replica boundaries.
    """
    frames = sorted(df['frame'].unique())
    unique_per_frame = []
    type_counts_per_frame = {itype: [] for itype in df['type'].unique()}

    for frame in frames:
        fd = df[df['frame'] == frame]
        unique_per_frame.append(fd['interaction_id'].nunique())
        for itype in type_counts_per_frame:
            type_counts_per_frame[itype].append(len(fd[fd['type'] == itype]))

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 10))

    ax1.plot(frames, unique_per_frame, marker='o', linewidth=2, markersize=4,
             color='#2E86AB', label='Unique Interactions')
    ax1.fill_between(frames, unique_per_frame, alpha=0.3, color='#2E86AB')

    window = max(3, len(frames) // 20)
    if len(unique_per_frame) >= window:
        rolling_avg = pd.Series(unique_per_frame).rolling(window=window, center=True).mean()
        ax1.plot(frames, rolling_avg, linewidth=2.5, color='#A23B72',
                 linestyle='--', label=f'Rolling Avg (window={window})')

    mean_val = np.mean(unique_per_frame)
    ax1.axhline(y=mean_val, color='red', linestyle=':', linewidth=1.5,
                label=f'Mean: {mean_val:.1f}')

    # Replica boundaries
    max_frame = frames[-1]
    boundaries = []
    b = replica_size
    while b < max_frame:
        boundaries.append(b)
        b += replica_size
    for bnd in boundaries:
        ax1.axvline(x=bnd, color='black', linestyle='--', linewidth=1.5, alpha=0.6)

    starts = [frames[0]] + boundaries
    ends   = boundaries + [max_frame]
    y_top  = max(25, max(unique_per_frame) + 2)
    for i, (s, e) in enumerate(zip(starts, ends)):
        ax1.text((s+e)/2, y_top*0.97, f'Rep{i+1}',
                 ha='center', va='top', fontsize=10, fontweight='bold',
                 color='darkgray', zorder=10)

    ax1.set_ylabel('Number of Unique Interactions', fontweight='bold', fontsize=12)
    ax1.set_title('Total Unique Interactions Over Time', fontweight='bold', fontsize=14, pad=15)
    ax1.legend(loc='lower left', framealpha=0.9)
    ax1.grid(True, alpha=0.3, linestyle='--')
    ax1.set_xlim(frames[0], frames[-1])
    ax1.set_ylim(0, y_top)

    type_colors = {
        'Hydrogen Bond': '#1f77b4', 'Hydrophobic': '#ff7f0e',
        'Salt Bridge': '#d62728',   'Pi-Pi Stacking': '#9467bd',
        'Pi-Cation': '#2ca02c',     'Halogen Bond': '#8c564b',
        'Metal Complex': '#e377c2'
    }
    stacked_data   = []
    stacked_labels = []
    stacked_colors = []
    for itype, counts in type_counts_per_frame.items():
        if sum(counts) > 0:
            stacked_data.append(counts)
            stacked_labels.append(itype)
            stacked_colors.append(type_colors.get(itype, 'gray'))

    ax2.stackplot(frames, *stacked_data, labels=stacked_labels,
                  colors=stacked_colors, alpha=0.8)
    ax2.set_xlabel('Frame Number', fontweight='bold', fontsize=12)
    ax2.set_ylabel('Number of Interactions', fontweight='bold', fontsize=12)
    ax2.set_title('Interaction Types Over Time (Stacked)', fontweight='bold', fontsize=14, pad=15)
    ax2.legend(loc='upper left', framealpha=0.9, fontsize=9)
    ax2.grid(True, alpha=0.3, linestyle='--')
    ax2.set_xlim(frames[0], frames[-1])

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_file}")

    print(f"\nUnique Interactions Statistics:")
    print(f"  Mean:    {np.mean(unique_per_frame):.1f}")
    print(f"  Std Dev: {np.std(unique_per_frame):.1f}")
    print(f"  Min:     {np.min(unique_per_frame)} (Frame {frames[np.argmin(unique_per_frame)]})")
    print(f"  Max:     {np.max(unique_per_frame)} (Frame {frames[np.argmax(unique_per_frame)]})")


# ============================================================
# ARGUMENT PARSER
# ============================================================

def parse_arguments():
    parser = argparse.ArgumentParser(
        description='PLIP Trajectory Analysis — Unified (manual + auto-detect)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
EXAMPLES:
  # Config file (set auto_detect_chains in config):
  python plip_trajectory_analysis.py config.yaml

  # CLI auto-detect mode:
  python plip_trajectory_analysis.py \\
      --topology system.prmtop \\
      --trajectories traj1.nc traj2.nc traj3.nc \\
      --output results/ --interval 32

  # CLI with replica size:
  python plip_trajectory_analysis.py \\
      --topology system.prmtop \\
      --trajectories rep1.nc rep2.nc rep3.nc \\
      --output results/ --replica-size 175
        """
    )
    parser.add_argument('config', nargs='?',
                        help='YAML config file (manual or auto-detect mode)')
    parser.add_argument('--topology',     help='Topology file (prmtop, pdb, etc)')
    parser.add_argument('--trajectories', nargs='+', help='Trajectory file(s)')
    parser.add_argument('--output',   default='plip_results',  help='Output directory')
    parser.add_argument('--prefix',   default='analysis',      help='Output file prefix')
    parser.add_argument('--start-frame',  type=int, default=1,   help='First frame (default: 1)')
    parser.add_argument('--end-frame',    type=int, default=-1,  help='Last frame (default: last)')
    parser.add_argument('--interval',     type=int, default=10,  help='Frame interval (default: 10)')
    parser.add_argument('--target-chain', default='B',
                        help='Chain to analyze interactions with (default: B)')
    parser.add_argument('--replica-size', type=int, default=175,
                        help='Frames per replica for timeline plot (default: 175)')
    return parser.parse_args()


# ============================================================
# MAIN PIPELINE
# ============================================================

def main():
    args = parse_arguments()

    # --- Route: config file or CLI ---
    if args.config and not args.topology:
        if not Path(args.config).exists():
            print(f"ERROR: Config file not found: {args.config}")
            sys.exit(1)
        config = Config(config_file=args.config)
    elif args.topology and args.trajectories:
        config = Config(args=args)
    else:
        print("ERROR: Provide either a config file OR --topology + --trajectories")
        print("Run with --help for usage information")
        sys.exit(1)

    # --- Header ---
    print(f"\n{'='*60}")
    print("PLIP TRAJECTORY ANALYSIS")
    print(f"{'='*60}")
    print(f"Topology:         {config.topology}")
    print(f"Trajectories ({len(config.trajectories)}):")
    for i, traj in enumerate(config.trajectories, 1):
        print(f"  {i}. {traj}")
    print(f"Frame interval:   {config.interval}")
    print(f"Output directory: {config.output_dir}")

    # --- Setup ---
    config.output_dir.mkdir(exist_ok=True, parents=True)
    print("\nChecking dependencies...")
    check_pdb_tools()
    if not shutil.which("cpptraj"):
        raise RuntimeError("cpptraj not found in PATH")
    if not shutil.which("plip"):
        raise RuntimeError("PLIP not found in PATH")

    work_dir = config.output_dir / "work"
    work_dir.mkdir(exist_ok=True)

    # --- Resolve chain assignment mode ---
    config.resolve_chain_mode(work_dir)

    # --- Extract frames ---
    frames = extract_frames(config, work_dir)
    total_frames = len(frames)

    # --- Process frames ---
    print(f"\n{'='*60}")
    print(f"PROCESSING {total_frames} FRAMES")
    print(f"{'='*60}")

    all_frame_data = {}
    failed_frames  = []
    skipped = new = 0

    for i, frame_pdb in enumerate(tqdm(frames, desc="Processing frames")):
        frame_num = i + 1
        processed_pdb  = work_dir / f"processed_frame_{frame_num}.pdb"
        plip_out_dir   = work_dir / f"plip_frame_{frame_num}"

        # Checkpoint: reuse existing PLIP XML if valid
        existing_xml = (list(plip_out_dir.glob("**/*report.xml"))
                        if plip_out_dir.exists() else [])
        if existing_xml:
            try:
                all_frame_data[frame_num] = parse_plip_xml(existing_xml[0])
                skipped += 1
                continue
            except Exception as e:
                tqdm.write(f"\n[CHECKPOINT] Frame {frame_num}: cached XML invalid ({e}), reprocessing...")

        try:
            fix_chain_assignments(frame_pdb, config)
            process_pdb_file(frame_pdb, processed_pdb, config)
            xml_file = run_plip(processed_pdb, config.target_chain, plip_out_dir)
            all_frame_data[frame_num] = parse_plip_xml(xml_file)
            new += 1
        except Exception as e:
            tqdm.write(f"\nWarning: Frame {frame_num} failed: {e}")
            failed_frames.append(frame_num)

    print(f"\n[CHECKPOINT SUMMARY]")
    print(f"  Skipped (cached): {skipped}")
    print(f"  Newly processed:  {new}")
    print(f"  Failed:           {len(failed_frames)}")
    print(f"  Total loaded:     {len(all_frame_data)}")

    if not all_frame_data:
        print("No frames were successfully processed!")
        sys.exit(1)
    if failed_frames:
        print(f"  Failed frames: {failed_frames}")

    # --- Aggregate ---
    print(f"\n{'='*60}")
    print("AGGREGATING INTERACTION DATA")
    print(f"{'='*60}")
    df = aggregate_interactions(all_frame_data)
    if df is None or df.empty:
        print("No interactions found across all frames!")
        sys.exit(1)
    prevalence_df = calculate_prevalence(df, len(all_frame_data))

    csv_file = config.output_dir / f"{config.output_prefix}_prevalence.csv"
    prevalence_df.to_csv(csv_file, index=False)
    print(f"\nSaved: {csv_file}")

    # --- Plots ---
    print(f"\n{'='*60}")
    print("GENERATING PLOTS")
    print(f"{'='*60}")
    pfx = config.output_dir / config.output_prefix

    plot_prevalence_overview(prevalence_df, f"{pfx}_prevalence_overview.png")
    plot_interaction_timeline(df, prevalence_df, f"{pfx}_timeline.png")
    plot_residue_heatmap(prevalence_df, f"{pfx}_residue_heatmap.png")
    plot_interaction_map(prevalence_df, f"{pfx}_interaction_map.png", min_prevalence=10)
    plot_chain_interaction_matrix(prevalence_df, f"{pfx}_chain_matrix.png")
    plot_interaction_fingerprint(df, prevalence_df, f"{pfx}_fingerprint.png",
                                 top_n=min(15, len(prevalence_df)))
    plot_unique_interactions_over_time(df, f"{pfx}_unique_over_time.png",
                                       replica_size=config.replica_size)

    # --- Summary ---
    print(f"\n{'='*60}")
    print("TRAJECTORY ANALYSIS SUMMARY")
    print(f"{'='*60}")
    print(f"Total frames analyzed:            {len(all_frame_data)}")
    print(f"Total unique interactions:         {len(prevalence_df)}")
    print(f"Average interactions per frame:    {len(df) / len(all_frame_data):.1f}")
    print(f"\nTop 10 most prevalent interactions:")
    for _, row in prevalence_df.head(10).iterrows():
        extra = f" [~{row['avg_per_frame']:.1f}x per frame]" if row.get('avg_per_frame', 1) > 1.1 else ""
        print(f"  {row['prot_res']}--{row['lig_res']} ({row['type']}): "
              f"{row['prevalence_pct']:.1f}% ({row['frame_count']}/{len(all_frame_data)} frames){extra}")

    print(f"\nAll results saved to: {config.output_dir}")
    print(f"\n{'='*60}")
    print("ANALYSIS COMPLETE!")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
