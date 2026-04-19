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

from Bio.PDB import PDBParser
from Bio.Align import PairwiseAligner
import subprocess

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

def _run_mkdssp_direct(pdb_path: str, mkdssp_path: str) -> str:
    """
    Call mkdssp directly via subprocess and return raw text output.
    Tries multiple command variants to handle mkdssp 4.x on Windows,
    where the mmcif_pdbx.dic dictionary may not be installed.
    """
    # Try these command variants in order until one produces output
    cmd_variants = [
        # Variant 1: explicit input + output format (best for mkdssp 4.x)
        [mkdssp_path, "--output-format", "dssp", "--input-format", "pdb", pdb_path],
        # Variant 2: output format only (mkdssp auto-detects input)
        [mkdssp_path, "--output-format", "dssp", pdb_path],
        # Variant 3: equals-sign syntax
        [mkdssp_path, "--output-format=dssp", pdb_path],
        # Variant 4: legacy syntax (mkdssp 3.x)
        [mkdssp_path, pdb_path],
    ]

    last_error = None
    for cmd in cmd_variants:
        print(f"[DSSP] Trying: {' '.join(cmd)}")
        result = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        if result.stderr.strip():
            print(f"[DSSP] stderr: {result.stderr.strip()[:200]}")

        if result.stdout.strip():
            print(f"[DSSP] Success with variant: {cmd[1]}")
            return result.stdout

        last_error = (cmd, result.returncode, result.stderr)

    # All variants failed — give actionable error message
    cmd, rc, stderr = last_error
    raise RuntimeError(
        f"mkdssp failed on all command variants.\n"
        f"Last command: {cmd}\n"
        f"Return code: {rc}\n"
        f"Stderr: {stderr[:500]}\n"
        f"\nFIX OPTIONS:\n"
        f"  1. Install missing dictionary: conda install -c conda-forge libcifpp\n"
        f"  2. Re-install DSSP:           conda install -c salilab dssp --force-reinstall\n"
        f"  3. Try a different conda env: conda create -n nmr python=3.10 && conda activate nmr\n"
        f"     then: conda install -c salilab dssp && pip install biopython pynmrstar pandas"
    )


def _parse_dssp_output(dssp_text: str, chain_id: Optional[str] = None) -> pd.DataFrame:
    """
    Parse raw mkdssp text output into a DataFrame.
    Handles the fixed-width DSSP format directly.
    """
    rows = []
    in_data = False

    for line in dssp_text.splitlines():
        if '#  RESIDUE' in line:
            in_data = True
            continue
        if not in_data or len(line) < 17:
            continue

        try:
            res_seq = line[0:5].strip()
            chain   = line[11].strip() if len(line) > 11 else ''
            aa      = line[13].strip() if len(line) > 13 else ''
            ss_raw  = line[16] if len(line) > 16 else ' '

            if not res_seq or aa in ('!', '*', ''):
                continue  # chain break / missing residue markers

            if chain_id and chain != chain_id:
                continue

            try:
                acc = int(line[34:38].strip()) if len(line) > 38 else None
            except ValueError:
                acc = None

            try:
                phi = float(line[103:109].strip()) if len(line) > 109 else None
                psi = float(line[109:115].strip()) if len(line) > 115 else None
            except ValueError:
                phi = psi = None

            rows.append({
                'residue_number': int(res_seq),
                'chain_id':       chain,
                'aa':             aa,
                'ss_raw':         ss_raw,
                'ss_class':       DSSP_RAW_TO_CLASS.get(ss_raw, 'coil'),
                'accessibility':  acc,
                'phi':            phi,
                'psi':            psi,
            })
        except (ValueError, IndexError):
            continue

    if not rows:
        raise RuntimeError("DSSP output parsed but no residue records found. Check PDB file.")

    return pd.DataFrame(rows).sort_values('residue_number').reset_index(drop=True)


def extract_dssp(pdb_path: str, chain_id: Optional[str] = None) -> Dict[int, str]:
    """
    Run DSSP on a PDB file.
    Returns: {residue_number (int) -> coarse_ss_class ('helix'/'strand'/'coil')}
    Uses direct subprocess call to avoid Biopython Windows stderr bug.
    """
    mkdssp_path = find_mkdssp()
    if mkdssp_path is None:
        raise EnvironmentError(
            "mkdssp executable not found.\n"
            "Fix: conda install -c salilab dssp"
        )

    dssp_text = _run_mkdssp_direct(pdb_path, mkdssp_path)
    df = _parse_dssp_output(dssp_text, chain_id=chain_id)
    ss_map = dict(zip(df['residue_number'], df['ss_class']))
    print(f"[DSSP] Extracted {len(ss_map)} residue assignments.")
    return ss_map


def extract_dssp_full(pdb_path: str, chain_id: Optional[str] = None) -> pd.DataFrame:
    """
    Run DSSP and return full DataFrame with columns:
        residue_number | chain_id | aa | ss_raw | ss_class | phi | psi | accessibility
    Uses direct subprocess call to avoid Biopython Windows stderr bug.
    """
    mkdssp_path = find_mkdssp()
    if mkdssp_path is None:
        raise EnvironmentError("mkdssp not found. Run: conda install -c salilab dssp")

    dssp_text = _run_mkdssp_direct(pdb_path, mkdssp_path)
    df = _parse_dssp_output(dssp_text, chain_id=chain_id)
    print(f"[DSSP] Extracted {len(df)} residues across chains: {sorted(df['chain_id'].unique())}")
    return df


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
# Pure-Python PDB fallback (when mkdssp is unavailable/broken)
# ---------------------------------------------------------------------------

def extract_ss_from_pdb_records(pdb_path: str) -> pd.DataFrame:
    """
    Extract secondary structure from HELIX and SHEET records in the PDB file.
    This is a fallback when mkdssp cannot run.

    Less detailed than DSSP (no phi/psi/accessibility) but works with zero
    external dependencies. PDB HELIX/SHEET records are deposited by the
    structure authors and are generally reliable.

    Returns DataFrame with same columns as extract_dssp_full.
    """
    helix_ranges = []   # list of (chain, start_res, end_res)
    sheet_ranges = []

    with open(pdb_path, 'r') as f:
        for line in f:
            rec = line[:6].strip()
            if rec == 'HELIX':
                try:
                    chain = line[19].strip()
                    start = int(line[21:25].strip())
                    end   = int(line[33:37].strip())
                    helix_ranges.append((chain, start, end))
                except (ValueError, IndexError):
                    continue
            elif rec == 'SHEET':
                try:
                    chain = line[21].strip()
                    start = int(line[22:26].strip())
                    end   = int(line[33:37].strip())
                    sheet_ranges.append((chain, start, end))
                except (ValueError, IndexError):
                    continue

    print(f"[PDB records] Found {len(helix_ranges)} helix ranges, {len(sheet_ranges)} sheet ranges")

    # Parse all residues from ATOM records
    residues = {}  # (chain, resnum) -> aa
    with open(pdb_path, 'r') as f:
        for line in f:
            if line[:4] != 'ATOM':
                continue
            try:
                chain  = line[21].strip()
                resnum = int(line[22:26].strip())
                resname = line[17:20].strip()
                key = (chain, resnum)
                if key not in residues:
                    THREE_TO_ONE = {
                        'ALA':'A','CYS':'C','ASP':'D','GLU':'E','PHE':'F','GLY':'G',
                        'HIS':'H','ILE':'I','LYS':'K','LEU':'L','MET':'M','ASN':'N',
                        'PRO':'P','GLN':'Q','ARG':'R','SER':'S','THR':'T','VAL':'V',
                        'TRP':'W','TYR':'Y',
                    }
                    residues[key] = THREE_TO_ONE.get(resname, 'X')
            except (ValueError, IndexError):
                continue

    if not residues:
        raise RuntimeError(f"No ATOM records found in {pdb_path}")

    def get_ss(chain, resnum):
        for (c, s, e) in helix_ranges:
            if c == chain and s <= resnum <= e:
                return 'H', 'helix'
        for (c, s, e) in sheet_ranges:
            if c == chain and s <= resnum <= e:
                return 'E', 'strand'
        return 'C', 'coil'

    rows = []
    for i, ((chain, resnum), aa) in enumerate(sorted(residues.items())):
        ss_raw, ss_class = get_ss(chain, resnum)
        rows.append({
            'residue_number': i + 1,   # sequential numbering like DSSP
            'pdb_resnum':     resnum,
            'chain_id':       chain,
            'aa':             aa,
            'ss_raw':         ss_raw,
            'ss_class':       ss_class,
            'accessibility':  None,
            'phi':            None,
            'psi':            None,
        })

    df = pd.DataFrame(rows)
    print(f"[PDB records] SS distribution: {df['ss_class'].value_counts().to_dict()}")
    return df


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
