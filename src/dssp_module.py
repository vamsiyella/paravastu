"""
dssp_module.py — PDB download, DSSP extraction, sequence alignment.

Windows-compatible: explicitly locates mkdssp in conda environment.
"""

import os
import sys
import shutil
import requests
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Optional, Dict, List, Tuple

from Bio.PDB import PDBParser, DSSP
from Bio.Align import PairwiseAligner

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

RCSB_URL = "https://files.rcsb.org/download/{pdb_id}.pdb"
EBI_URL  = "https://www.ebi.ac.uk/pdbe/entry-files/download/pdb{pdb_id_lower}.ent"

DSSP_RAW_TO_CLASS = {
    'H': 'helix',   # alpha helix
    'G': 'helix',   # 3-10 helix
    'I': 'helix',   # pi helix
    'E': 'strand',  # beta strand
    'B': 'strand',  # beta bridge
    'T': 'coil',    # turn
    'S': 'coil',    # bend
    'C': 'coil',    # coil / unassigned
    '-': 'coil',    # sometimes used for unassigned
    ' ': 'coil',
}


# ---------------------------------------------------------------------------
# DSSP executable detection (critical for Windows)
# ---------------------------------------------------------------------------

def find_mkdssp() -> Optional[str]:
    """
    Locate mkdssp executable.
    Search order:
      1. PATH (works if conda env is active)
      2. Conda envs directory (handles Windows where PATH may not be set)
      3. Common install locations
    Returns absolute path string, or None if not found.
    """
    # 1. Check PATH first
    found = shutil.which("mkdssp") or shutil.which("dssp")
    if found:
        return found

    # 2. Search conda environments
    conda_root = os.environ.get("CONDA_PREFIX") or os.environ.get("CONDA_ROOT")
    if conda_root:
        candidates = [
            Path(conda_root) / "bin" / "mkdssp",
            Path(conda_root) / "Library" / "bin" / "mkdssp.exe",  # Windows conda
            Path(conda_root) / "Scripts" / "mkdssp.exe",
        ]
        for c in candidates:
            if c.exists():
                return str(c)

    # 3. Walk conda envs if CONDA_PREFIX not set
    home = Path.home()
    possible_roots = [
        home / "anaconda3",
        home / "miniconda3",
        home / "anaconda",
        home / "miniconda",
        Path("C:/ProgramData/Anaconda3"),
        Path("C:/Users") / os.environ.get("USERNAME", "") / "anaconda3",
        Path("C:/Users") / os.environ.get("USERNAME", "") / "miniconda3",
    ]
    for root in possible_roots:
        for pattern in ["**/mkdssp", "**/mkdssp.exe"]:
            matches = list(root.glob(pattern))
            if matches:
                return str(matches[0])

    return None


# ---------------------------------------------------------------------------
# PDB download
# ---------------------------------------------------------------------------

def download_pdb(pdb_id: str, data_dir: Path) -> str:
    """
    Download PDB file from RCSB (or EBI as fallback).
    Returns local file path string.
    """
    data_dir = Path(data_dir)
    data_dir.mkdir(parents=True, exist_ok=True)

    pdb_path = data_dir / f"{pdb_id}.pdb"
    if pdb_path.exists():
        print(f"[PDB] Using cached: {pdb_path}")
        return str(pdb_path)

    urls = [
        RCSB_URL.format(pdb_id=pdb_id),
        EBI_URL.format(pdb_id_lower=pdb_id.lower()),
    ]

    for url in urls:
        try:
            print(f"[PDB] Downloading {url}")
            resp = requests.get(url, timeout=30)
            resp.raise_for_status()
            pdb_path.write_text(resp.text)
            print(f"[PDB] Saved to {pdb_path}")
            return str(pdb_path)
        except Exception as e:
            print(f"[PDB] Failed {url}: {e}")
            continue

    raise RuntimeError(
        f"Could not download PDB {pdb_id}. "
        "Place {pdb_id}.pdb manually in the data/ directory and retry."
    )


# ---------------------------------------------------------------------------
# DSSP extraction
# ---------------------------------------------------------------------------

def extract_dssp(pdb_path: str, chain_id: Optional[str] = None) -> Dict[int, str]:
    """
    Run DSSP on a PDB file.
    Returns: {residue_number (int) -> coarse_ss_class ('helix'/'strand'/'coil')}

    chain_id: restrict to this chain. If None, uses first chain found.
    """
    mkdssp_path = find_mkdssp()
    if mkdssp_path is None:
        raise EnvironmentError(
            "mkdssp executable not found.\n"
            "Fix: conda install -c salilab dssp\n"
            "Or install from: https://swift.cmbi.umcn.nl/gv/dssp/"
        )
    print(f"[DSSP] Using executable: {mkdssp_path}")

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_path)
    model = structure[0]

    if chain_id is None:
        chains = list(model.get_chains())
        if not chains:
            raise ValueError("No chains found in PDB structure.")
        chain_id = chains[0].id
        print(f"[DSSP] Using chain: {chain_id}")

    # Run DSSP - pass executable path explicitly for Windows compatibility
    try:
        dssp = DSSP(model, pdb_path, dssp=mkdssp_path)
    except Exception as e:
        raise RuntimeError(f"DSSP failed: {e}\nPDB: {pdb_path}\nExecutable: {mkdssp_path}")

    ss_map = {}
    for key in dssp.keys():
        ch, res_id = key
        if ch != chain_id:
            continue
        residue_number = res_id[1]
        ss_raw = dssp[key][2]  # character code
        ss_class = DSSP_RAW_TO_CLASS.get(ss_raw, 'coil')
        ss_map[residue_number] = ss_class

    print(f"[DSSP] Extracted {len(ss_map)} residue assignments.")
    return ss_map


def extract_dssp_full(pdb_path: str, chain_id: Optional[str] = None) -> pd.DataFrame:
    """
    Run DSSP and return full DataFrame with columns:
        residue_number | chain_id | ss_raw | ss_class
    """
    mkdssp_path = find_mkdssp()
    if mkdssp_path is None:
        raise EnvironmentError("mkdssp not found. Run: conda install -c salilab dssp")

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_path)
    model = structure[0]

    if chain_id is None:
        chain_id = list(model.get_chains())[0].id

    dssp = DSSP(model, pdb_path, dssp=mkdssp_path)

    rows = []
    for key in dssp.keys():
        ch, res_id = key
        if ch != chain_id:
            continue
        resnum = res_id[1]
        dssp_data = dssp[key]
        ss_raw   = dssp_data[2]
        aa       = dssp_data[1]
        phi      = dssp_data[4]
        psi      = dssp_data[5]
        acc      = dssp_data[3]  # solvent accessibility
        rows.append({
            'residue_number': resnum,
            'chain_id':       ch,
            'aa':             aa,
            'ss_raw':         ss_raw,
            'ss_class':       DSSP_RAW_TO_CLASS.get(ss_raw, 'coil'),
            'phi':            phi,
            'psi':            psi,
            'accessibility':  acc,
        })

    return pd.DataFrame(rows).sort_values('residue_number').reset_index(drop=True)


# ---------------------------------------------------------------------------
# Secondary structure segmentation
# ---------------------------------------------------------------------------

def build_ss_segments(ss_map: Dict[int, str]) -> List[dict]:
    """
    Convert per-residue SS map into contiguous segments.
    Returns list of {start, end, ss_class, length}.
    """
    segments = []
    sorted_res = sorted(ss_map.keys())
    if not sorted_res:
        return segments

    start = sorted_res[0]
    current_ss = ss_map[start]
    prev = start

    for res in sorted_res[1:]:
        if res == prev + 1 and ss_map[res] == current_ss:
            prev = res
            continue
        segments.append({'start': start, 'end': prev, 'ss_class': current_ss, 'length': prev - start + 1})
        start = res
        current_ss = ss_map[res]
        prev = res

    segments.append({'start': start, 'end': prev, 'ss_class': current_ss, 'length': prev - start + 1})
    return segments


# ---------------------------------------------------------------------------
# Sequence alignment for BMRB <-> PDB residue number mapping
# ---------------------------------------------------------------------------

def align_and_map(bmrb_seq: str, pdb_seq: str) -> Dict[int, int]:
    """
    Align BMRB sequence to PDB sequence using global pairwise alignment.
    Returns mapping: bmrb_seq_id (1-based) -> pdb_residue_number (1-based).

    This handles cases where:
    - Signal peptides are cleaved in the PDB structure
    - BMRB and PDB have different numbering offsets
    - Partial sequences
    """
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -0.5

    alignments = list(aligner.align(bmrb_seq, pdb_seq))

    if not alignments:
        print("[ALIGN] No alignment found. Using identity mapping.")
        return {i+1: i+1 for i in range(len(bmrb_seq))}

    best = alignments[0]
    aligned_bmrb = best[0]
    aligned_pdb  = best[1]

    # Build position maps (skip gaps)
    bmrb_positions = []
    pos = 0
    for ch in aligned_bmrb:
        if ch != '-':
            pos += 1
            bmrb_positions.append(pos)
        else:
            bmrb_positions.append(None)

    pdb_positions = []
    pos = 0
    for ch in aligned_pdb:
        if ch != '-':
            pos += 1
            pdb_positions.append(pos)
        else:
            pdb_positions.append(None)

    mapping = {}
    for bp, pp in zip(bmrb_positions, pdb_positions):
        if bp is not None and pp is not None:
            mapping[bp] = pp

    print(f"[ALIGN] Mapped {len(mapping)} residues (BMRB seq len {len(bmrb_seq)}, PDB seq len {len(pdb_seq)})")
    return mapping


def get_pdb_sequence_from_file(pdb_path: str, chain_id: Optional[str] = None) -> str:
    """Extract amino acid sequence from a PDB file."""
    THREE_TO_ONE = {
        'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
        'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
        'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
        'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y',
    }
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_path)
    model = structure[0]

    if chain_id is None:
        chain_id = list(model.get_chains())[0].id

    chain = model[chain_id]
    seq = ''
    for residue in chain.get_residues():
        resname = residue.get_resname().strip()
        if resname in THREE_TO_ONE:
            seq += THREE_TO_ONE[resname]
    return seq
