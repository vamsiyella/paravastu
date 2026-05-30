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
    Returns local file path string (prefers .cif for mkdssp 4.x compatibility).
    """
    data_dir = Path(data_dir)
    data_dir.mkdir(parents=True, exist_ok=True)

    # Prefer .cif — mkdssp 4.x requires mmCIF format (PDB triggers dictionary error)
    cif_path = data_dir / f"{pdb_id}.cif"
    pdb_path = data_dir / f"{pdb_id}.pdb"

    if cif_path.exists():
        print(f"[PDB] Using cached CIF: {cif_path}")
        return str(cif_path)
    if pdb_path.exists():
        print(f"[PDB] Using cached PDB: {pdb_path}")
        return str(pdb_path)

    # Try CIF first (native mkdssp 4.x format)
    cif_urls = [
        f"https://files.rcsb.org/download/{pdb_id}.cif",
        f"https://www.ebi.ac.uk/pdbe/entry-files/download/{pdb_id.lower()}.cif",
    ]
    for url in cif_urls:
        try:
            print(f"[PDB] Downloading CIF: {url}")
            resp = requests.get(url, timeout=30)
            resp.raise_for_status()
            cif_path.write_text(resp.text)
            print(f"[PDB] Saved CIF to {cif_path}")
            return str(cif_path)
        except Exception as e:
            print(f"[PDB] CIF failed: {e}")

    # Fall back to PDB format
    pdb_urls = [
        RCSB_URL.format(pdb_id=pdb_id),
        EBI_URL.format(pdb_id_lower=pdb_id.lower()),
    ]
    for url in pdb_urls:
        try:
            print(f"[PDB] Downloading PDB: {url}")
            resp = requests.get(url, timeout=30)
            resp.raise_for_status()
            pdb_path.write_text(resp.text)
            print(f"[PDB] Saved PDB to {pdb_path}")
            return str(pdb_path)
        except Exception as e:
            print(f"[PDB] PDB failed: {e}")

    raise RuntimeError(
        f"Could not download {pdb_id} in any format.\n"
        f"Manually download from https://www.rcsb.org/structure/{pdb_id}\n"
        f"Save as data/{pdb_id}.cif (preferred) or data/{pdb_id}.pdb"
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
    # Detect file format from extension
    is_cif = str(pdb_path).lower().endswith('.cif')

    # Try these command variants in order until one produces output
    if is_cif:
        # mmCIF: mkdssp 4.x native format — no --input-format flag needed
        cmd_variants = [
            [mkdssp_path, "--output-format", "dssp", pdb_path],
            [mkdssp_path, "--output-format=dssp", pdb_path],
            [mkdssp_path, pdb_path],
        ]
    else:
        # PDB format: try explicit flag first, then fallbacks
        cmd_variants = [
            [mkdssp_path, "--output-format", "dssp", "--input-format", "pdb", pdb_path],
            [mkdssp_path, "--output-format", "dssp", pdb_path],
            [mkdssp_path, "--output-format=dssp", pdb_path],
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
    Extract secondary structure from HELIX/SHEET records (PDB) or
    _struct_conf/_struct_sheet_range (mmCIF).
    Fallback when mkdssp cannot run.

    Returns DataFrame with same columns as extract_dssp_full.
    """
    helix_ranges = []   # list of (chain, start_res, end_res)
    sheet_ranges = []
    is_cif = str(pdb_path).lower().endswith('.cif')

    if is_cif:
        # --- mmCIF: parse _struct_conf (helices) and _struct_sheet_range (strands) ---
        current_category = None
        col_names = []
        in_loop = False

        with open(pdb_path, 'r') as f:
            for line in f:
                line = line.rstrip()

                if line.startswith('loop_'):
                    in_loop = True
                    col_names = []
                    current_category = None
                    continue

                if line.startswith('_struct_conf.'):
                    current_category = 'struct_conf'
                    col = line.split('.')[1].strip().split()[0]
                    col_names.append(col)
                    continue

                if line.startswith('_struct_sheet_range.'):
                    current_category = 'struct_sheet_range'
                    col = line.split('.')[1].strip().split()[0]
                    col_names.append(col)
                    continue

                if line.startswith('_') and current_category:
                    # Different category started — reset
                    current_category = None
                    col_names = []
                    in_loop = False
                    continue

                if not line or line.startswith('#'):
                    continue

                if current_category and col_names and in_loop:
                    parts = line.split()
                    if not parts or parts[0].startswith('_'):
                        continue
                    ci = {c: i for i, c in enumerate(col_names)}

                    if current_category == 'struct_conf':
                        try:
                            conf_type = parts[ci['conf_type_id']]
                            chain = parts[ci['beg_label_asym_id']]
                            start = end = None
                            for sk, ek in [('beg_label_seq_id','end_label_seq_id'),
                                        ('beg_auth_seq_id','end_auth_seq_id')]:
                                if sk in ci and ek in ci:
                                    try:
                                        start = int(parts[ci[sk]])
                                        end   = int(parts[ci[ek]])
                                        break
                                    except (ValueError, IndexError):
                                        pass
                            if start is None:
                                continue
                            if 'HELX' in conf_type.upper():
                                helix_ranges.append((chain, start, end))
                            elif 'STRN' in conf_type.upper():   # beta-solenoid strands
                                sheet_ranges.append((chain, start, end))
                        except (KeyError, ValueError, IndexError):
                            pass

                    elif current_category == 'struct_sheet_range':
                        try:
                            chain = parts[ci['beg_label_asym_id']]
                            start = end = None
                            for sk, ek in [('beg_label_seq_id','end_label_seq_id'),
                                        ('beg_auth_seq_id','end_auth_seq_id')]:
                                if sk in ci and ek in ci:
                                    try:
                                        start = int(parts[ci[sk]])
                                        end   = int(parts[ci[ek]])
                                        break
                                    except (ValueError, IndexError):
                                        pass
                            if start is None:
                                continue
                            sheet_ranges.append((chain, start, end))
                        except (KeyError, ValueError, IndexError):
                            pass
    else:
        # --- PDB format: HELIX and SHEET records ---
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
    # If CIF gave nothing, try the .pdb counterpart (some NMR structures
    # have HELIX/SHEET in .pdb but not in .cif)
    if len(helix_ranges) + len(sheet_ranges) == 0:
        pdb_counterpart = str(pdb_path).replace('.cif', '.pdb')
        if Path(pdb_counterpart).exists():
            print(f"[PDB records] CIF empty → trying .pdb: {pdb_counterpart}")
            with open(pdb_counterpart, 'r') as f:
                for line in f:
                    rec = line[:6].strip()
                    if rec == 'HELIX':
                        try:
                            helix_ranges.append((line[19].strip(),
                                int(line[21:25]), int(line[33:37])))
                        except (ValueError, IndexError): pass
                    elif rec == 'SHEET':
                        try:
                            sheet_ranges.append((line[21].strip(),
                                int(line[22:26]), int(line[33:37])))
                        except (ValueError, IndexError): pass
            print(f"[PDB records] .pdb fallback: {len(helix_ranges)} helix, {len(sheet_ranges)} sheet")
            
    THREE_TO_ONE = {
        'ALA':'A','CYS':'C','ASP':'D','GLU':'E','PHE':'F','GLY':'G',
        'HIS':'H','ILE':'I','LYS':'K','LEU':'L','MET':'M','ASN':'N',
        'PRO':'P','GLN':'Q','ARG':'R','SER':'S','THR':'T','VAL':'V',
        'TRP':'W','TYR':'Y',
    }

    # Parse all residues from ATOM records
    residues = {}  # (chain, seq_id_int) -> aa

    if is_cif:
        # Parse _atom_site loop directly — more reliable than Biopython for minimal CIFs
        atom_cols = []
        in_atom_loop = False
        with open(pdb_path, 'r') as f:
            for line in f:
                line = line.rstrip()
                if line.startswith('_atom_site.'):
                    col = line.split('.')[1].strip().split()[0]
                    atom_cols.append(col)
                    in_atom_loop = True
                    continue
                if in_atom_loop and not line.startswith('_') and not line.startswith('#') and line.strip():
                    parts = line.split()
                    if not parts or parts[0].startswith('_') or parts[0] == 'loop_':
                        in_atom_loop = False
                        continue
                    if not atom_cols:
                        continue
                    ci = {c: i for i, c in enumerate(atom_cols)}
                    try:
                        group   = parts[ci.get('group_PDB', 0)]
                        if group not in ('ATOM', 'ATOM?'):
                            continue
                        chain   = parts[ci['label_asym_id']]
                        seq_id  = int(parts[ci['label_seq_id']])
                        resname = parts[ci['label_comp_id']]
                        key = (chain, seq_id)
                        if key not in residues:
                            residues[key] = THREE_TO_ONE.get(resname, 'X')
                    except (KeyError, ValueError, IndexError):
                        continue
    else:
        # PDB format: fixed-width ATOM records
        with open(pdb_path, 'r') as f:
            for line in f:
                if line[:4] != 'ATOM':
                    continue
                try:
                    chain   = line[21].strip()
                    resnum  = int(line[22:26].strip())
                    resname = line[17:20].strip()
                    key = (chain, resnum)
                    if key not in residues:
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
# Build BMRB-indexed SS map for predictions
# ---------------------------------------------------------------------------

def build_prediction_ss_map(
    dssp_df: pd.DataFrame,
    seq_mapping: Dict[int, int],
    sequence_length: int,
) -> Dict[int, str]:
    """
    Build a {bmrb_seq_id -> ss_class} map for use with ShiftPredictor.

    dssp_df:        output of extract_dssp_full or extract_ss_from_pdb_records
                    must have columns: pdb_resnum (or residue_number), ss_class
    seq_mapping:    {bmrb_seq_id -> pdb_aligned_position (1-based)}
                    output of align_and_map()
    sequence_length: length of the BMRB sequence

    Without this, predictions use the raw DSSP sequential numbering which
    doesn't match BMRB residue numbering when PDB has multiple chains or
    non-sequential residue numbers.
    """
    # Build pdb_position -> ss_class lookup
    # Use pdb_resnum if available, else fall back to residue_number
    pos_col = 'pdb_resnum' if 'pdb_resnum' in dssp_df.columns else 'residue_number'

    # For multi-chain PDB (like 2LBH with 2 identical chains), restrict to chain A
    if 'chain_id' in dssp_df.columns and dssp_df['chain_id'].nunique() > 1:
        first_chain = sorted(dssp_df['chain_id'].unique())[0]
        df_chain = dssp_df[dssp_df['chain_id'] == first_chain].copy()
        print(f"[SS MAP] Multi-chain PDB — using chain {first_chain} for mapping")
    else:
        df_chain = dssp_df.copy()

    # Map: pdb_position -> ss_class
    # Re-index as sequential position within the chosen chain
    df_chain = df_chain.reset_index(drop=True)
    df_chain['chain_pos'] = range(1, len(df_chain) + 1)
    pos_to_ss = dict(zip(df_chain['chain_pos'], df_chain['ss_class']))

    # Apply seq_mapping: bmrb_id -> pdb_position -> ss_class
    ss_map = {}
    for bmrb_id in range(1, sequence_length + 1):
        pdb_pos = seq_mapping.get(bmrb_id)
        if pdb_pos is not None:
            ss_map[bmrb_id] = pos_to_ss.get(pdb_pos, 'coil')
        else:
            ss_map[bmrb_id] = 'coil'  # unmapped residues default to coil

    mapped = sum(1 for v in ss_map.values() if v != 'coil')
    ss_counts = {}
    for v in ss_map.values():
        ss_counts[v] = ss_counts.get(v, 0) + 1
    print(f"[SS MAP] Built prediction map: {ss_counts}")
    return ss_map


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
    """Extract amino acid sequence from a PDB or mmCIF file."""
    THREE_TO_ONE = {
        'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
        'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
        'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
        'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y',
    }

    is_cif = str(pdb_path).lower().endswith('.cif')

    if is_cif:
        from Bio.PDB.MMCIFParser import MMCIFParser
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)

    structure = parser.get_structure("protein", pdb_path)
    model = structure[0]

    if chain_id is None:
        chain_id = sorted(c.id for c in model.get_chains())[0]

    chain = model[chain_id]
    seq = ''
    for residue in chain.get_residues():
        resname = residue.get_resname().strip()
        if resname in THREE_TO_ONE:
            seq += THREE_TO_ONE[resname]
    return seq
