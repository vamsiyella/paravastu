"""
bmrb_module.py — Fetch and parse BMRB NMR-STAR entries.

Handles:
- Downloading NMR-STAR files from BMRB
- Parsing chemical shift loops
- Sequence extraction from shift assignments
- PDB ID extraction from metadata
"""

import re
import requests
import pynmrstar
import pandas as pd
from pathlib import Path
from collections import defaultdict
from typing import Optional

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

BMRB_URL_TEMPLATE = (
    "https://bmrb.io/ftp/pub/bmrb/entry_directories/bmr{id}/bmr{id}_3.str"
)
BMRB_ALT_URL_TEMPLATE = (
    "https://bmrb.io/ftp/pub/bmrb/entry_directories/bmr{id}/bmr{id}.str"
)

THREE_TO_ONE = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
    'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
    'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
    'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y',
}

ATOM_TYPES = {'C', 'CA', 'CB', 'N', 'H', 'HA', 'HB', 'HB2', 'HB3', 'CO'}


# ---------------------------------------------------------------------------
# Fetching
# ---------------------------------------------------------------------------

def fetch_bmrb_text(bmr_id: int, cache_dir: Optional[Path] = None) -> str:
    """
    Download NMR-STAR text for a BMRB entry.
    Caches to disk if cache_dir is provided.
    Tries both _3.str and .str URL variants.
    """
    if cache_dir is not None:
        cache_dir = Path(cache_dir)
        cache_dir.mkdir(parents=True, exist_ok=True)
        cache_path = cache_dir / f"bmr{bmr_id}_3.str"
        if cache_path.exists():
            print(f"[BMRB] Loading cached: {cache_path}")
            return cache_path.read_text()

    urls = [
        BMRB_URL_TEMPLATE.format(id=bmr_id),
        BMRB_ALT_URL_TEMPLATE.format(id=bmr_id),
    ]
    last_error = None
    for url in urls:
        try:
            print(f"[BMRB] Fetching {url}")
            resp = requests.get(url, timeout=30)
            resp.raise_for_status()
            text = resp.text
            if cache_dir is not None:
                cache_path.write_text(text)
                print(f"[BMRB] Cached to {cache_path}")
            return text
        except Exception as e:
            last_error = e
            continue

    raise RuntimeError(
        f"Failed to fetch BMRB {bmr_id}. Last error: {last_error}\n"
        "If running locally, check network. "
        "You can also place the .str file manually in data/ and use load_bmrb_from_file()."
    )


def load_bmrb_from_file(path: str | Path) -> str:
    """Load NMR-STAR text from a local file (fallback when network unavailable)."""
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"BMRB file not found: {path}")
    return path.read_text()


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

def parse_shifts(nmrstar_text: str) -> pd.DataFrame:
    """
    Parse Atom_chem_shift loop from NMR-STAR text.
    Returns DataFrame with columns:
        seq_id, residue_name, residue, atom, shift
    """
    entry = pynmrstar.Entry.from_string(nmrstar_text)
    loops = entry.get_loops_by_category('Atom_chem_shift')

    if not loops:
        print("[BMRB] Warning: no Atom_chem_shift loops found.")
        return pd.DataFrame()

    rows = []
    for loop in loops:
        tags = loop.tags
        required = {'Comp_ID', 'Atom_ID', 'Val', 'Seq_ID'}
        if not required.issubset(set(tags)):
            print(f"[BMRB] Skipping loop — missing required tags. Found: {tags}")
            continue

        idx = {tag: i for i, tag in enumerate(tags)}

        for row in loop.data:
            try:
                comp    = str(row[idx['Comp_ID']]).strip()
                atom    = str(row[idx['Atom_ID']]).strip()
                val_str = str(row[idx['Val']]).strip()
                seq_id  = str(row[idx['Seq_ID']]).strip()

                if val_str in ('.', '?', ''):
                    continue
                shift  = float(val_str)
                seq_id = int(seq_id)

            except (ValueError, IndexError, TypeError):
                continue

            rows.append({
                'seq_id':       seq_id,
                'residue_name': comp,
                'residue':      THREE_TO_ONE.get(comp.upper(), 'X'),
                'atom':         atom,
                'shift':        shift,
            })

    df = pd.DataFrame(rows)
    if not df.empty:
        df = df.sort_values(['seq_id', 'atom']).reset_index(drop=True)
    return df


def extract_sequence(nmrstar_text: str) -> Optional[str]:
    """
    Extract amino acid sequence from NMR-STAR.
    Strategy: infer from Atom_chem_shift loop (works even when no
    Entity_poly_seq is present, which is common for ssNMR entries).
    """
    entry = pynmrstar.Entry.from_string(nmrstar_text)
    loops = entry.get_loops_by_category('Atom_chem_shift')

    seq_records = {}
    for loop in loops:
        tags = loop.tags
        if 'Comp_ID' not in tags or 'Seq_ID' not in tags:
            continue
        comp_idx = tags.index('Comp_ID')
        seq_idx  = tags.index('Seq_ID')

        for row in loop.data:
            try:
                resnum  = int(row[seq_idx])
                resname = str(row[comp_idx]).strip()
            except (TypeError, ValueError, IndexError):
                continue
            if resnum not in seq_records:
                seq_records[resnum] = resname

    if not seq_records:
        return None

    return ''.join(
        THREE_TO_ONE.get(seq_records[i].upper(), 'X')
        for i in sorted(seq_records)
    )


def extract_pdb_id(nmrstar_text: str) -> Optional[str]:
    """
    Try to extract an associated PDB ID from BMRB metadata.
    Checks multiple saveframe categories.
    """
    entry = pynmrstar.Entry.from_string(nmrstar_text)

    pdb_tag_names = ['PDB_ID', 'Database_accession_code', 'Related_BMRB_accession_code']

    for sf in entry.frame_list:
        for tag in pdb_tag_names:
            try:
                values = sf.get_tag(tag)
                for v in values:
                    v = str(v).strip()
                    # PDB IDs are 4 chars, alphanumeric
                    if re.match(r'^[A-Za-z0-9]{4}$', v):
                        return v.upper()
            except (KeyError, AttributeError):
                continue

    return None


# ---------------------------------------------------------------------------
# Coverage analysis
# ---------------------------------------------------------------------------

def compute_coverage(df: pd.DataFrame) -> pd.DataFrame:
    """
    For each residue, report which atoms are observed and how many.
    Returns DataFrame: seq_id | residue | atoms_observed | n_atoms
    """
    grp = df.groupby(['seq_id', 'residue'])['atom'].apply(
        lambda x: sorted(set(x))
    ).reset_index()
    grp.columns = ['seq_id', 'residue', 'atoms_observed']
    grp['n_atoms'] = grp['atoms_observed'].apply(len)
    return grp.sort_values('seq_id').reset_index(drop=True)


def compute_shift_stats(df: pd.DataFrame) -> pd.DataFrame:
    """
    Aggregate shift statistics grouped by (residue, atom, ss_class).
    ss_class column must exist; use 'unknown' if DSSP not yet run.

    Returns DataFrame with:
        residue | atom | ss_class | count | mean | median | std | min | max
    """
    if 'ss_class' not in df.columns:
        df = df.copy()
        df['ss_class'] = 'unknown'

    grouped = df.groupby(['residue', 'atom', 'ss_class'])['shift']
    stats = grouped.agg(
        count='count',
        mean='mean',
        median='median',
        std='std',
        min='min',
        max='max',
    ).reset_index()
    return stats.sort_values(['residue', 'atom', 'ss_class']).reset_index(drop=True)
