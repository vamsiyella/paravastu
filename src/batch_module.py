"""
batch_module.py — Phase 3: Scale the Reference Database

Curated list of solid-state NMR BMRB entries with linked PDB structures.
All entries are beta-sheet-rich or mixed topology proteins characterized
by solid-state NMR — the exact use case for this pipeline.

Batch runner processes all entries and merges their shift statistics into
a single reference database split by helix/strand/coil.
"""

import pandas as pd
import numpy as np
from pathlib import Path
from collections import defaultdict
import json
import time

# ---------------------------------------------------------------------------
# Curated solid-state BMRB entry list
# Format: (bmrb_id, pdb_id, description, topology_notes)
#
# Selection criteria:
#   - Solid-state NMR (not solution NMR)
#   - 13C and/or 15N chemical shifts deposited
#   - Linked PDB structure available
#   - Well-characterized, published entries
# ---------------------------------------------------------------------------

SOLID_STATE_ENTRIES = [
    # ---- Amyloid fibrils ----
    (17661, "2LMP", "Ure2p prion domain fibril — mixed coil/strand"),
    (15409, "2KIB", "HET-s prion domain — beta-solenoid, all strand"),
    (18525, "2M5N", "alpha-synuclein fibril — strand-rich"),
    (30430, "6SST", "tau fibril — beta-sheet core"),
    (26612, "5KK3", "amyloid-beta 1-42 fibril — all strand"),

    # ---- Membrane proteins ----
    (6838,  "1YMZ", "fd coat protein — helix-rich membrane anchor"),
    (17557, "2KSF", "M2 proton channel — helix bundle"),
    (15523, "2L9B", "gramicidin A — beta-helix"),

    # ---- Cysteine-rich / structured peptides ----
    (17561, "2LBH", "EETI-II knottin — beta-sheet + coil (Phase 1/2 entry)"),
    (16318, "2JWU", "ubiquitin microcrystals — alpha/beta mixed"),
    (15380, "1LY2", "GB1 protein — alpha/beta mixed"),

    # ---- Beta-sheet rich proteins ----
    (15865, "1M8M", "SH3 domain — all beta"),
    (17948, "2NUZ", "calmodulin — helix-rich"),
    (5969,  "1H4W", "BPTI — disulfide-rich beta"),

    # ---- Mixed topology ----
    (16299, "2JSV", "thioredoxin microcrystals — alpha/beta"),
    (17604, "2KXA", "DsbA N-terminal domain — mixed"),
]

# Entries known to have incomplete or incompatible data — skip these
SKIP_ENTRIES = {
    30430,  # tau — PDB 6SST requires special chain handling
    26612,  # Abeta fibril — multiple incompatible models
}

# ---------------------------------------------------------------------------
# Batch runner
# ---------------------------------------------------------------------------

def run_batch(
    bmrb_pdb_pairs=None,
    output_dir=None,
    skip_ids=None,
    verbose=True,
    cache_dir=None,
):
    """
    Process multiple BMRB entries and aggregate shift statistics.

    Parameters
    ----------
    bmrb_pdb_pairs : list of (bmrb_id, pdb_id) or None
        If None, uses the built-in SOLID_STATE_ENTRIES list.
        Pass your own list to override, e.g. [(17561, '2LBH'), (15409, '2KIB')]
    output_dir : str or Path
        Where to save the merged reference CSV.
        Defaults to <repo_root>/data/
    skip_ids : set of int, optional
        BMRB IDs to skip (merged with SKIP_ENTRIES).
    verbose : bool
        Print per-entry progress.
    cache_dir : str or Path, optional
        Where raw per-entry CSVs are cached. Defaults to output_dir/batch_cache/

    Returns
    -------
    pd.DataFrame
        Merged statistics table: residue | atom | ss_class | count | mean | ...
        Saved to output_dir/reference_db.csv
    """
    from bmrb_module import fetch_bmrb_star, parse_shifts, extract_sequence_from_bmrb
    from dssp_module import download_cif, extract_ss_from_cif, download_pdb, align_sequences_and_map

    if bmrb_pdb_pairs is None:
        bmrb_pdb_pairs = [(e[0], e[1]) for e in SOLID_STATE_ENTRIES]

    if skip_ids is None:
        skip_ids = set()
    skip_ids = skip_ids | SKIP_ENTRIES

    # Resolve paths
    repo_root = Path(__file__).resolve().parent.parent
    if output_dir is None:
        output_dir = repo_root / "data"
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)

    if cache_dir is None:
        cache_dir = output_dir / "batch_cache"
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(exist_ok=True)

    # ---------- process each entry ----------
    all_raw_shifts = []   # list of DataFrames, one per entry
    results_log = []

    for bmrb_id, pdb_id in bmrb_pdb_pairs:
        if bmrb_id in skip_ids:
            if verbose:
                print(f"  [{bmrb_id}/{pdb_id}] SKIPPED (in skip list)")
            results_log.append({"bmrb": bmrb_id, "pdb": pdb_id, "status": "skipped", "shifts": 0})
            continue

        cache_file = cache_dir / f"bmr{bmrb_id}_{pdb_id}_raw.csv"

        if cache_file.exists():
            if verbose:
                print(f"  [{bmrb_id}/{pdb_id}] Loading from cache...")
            try:
                df = pd.read_csv(cache_file)
                all_raw_shifts.append(df)
                results_log.append({"bmrb": bmrb_id, "pdb": pdb_id, "status": "cached", "shifts": len(df)})
                continue
            except Exception as e:
                if verbose:
                    print(f"    Cache read failed ({e}), re-processing...")

        if verbose:
            print(f"\n{'='*60}")
            print(f"  Processing BMRB {bmrb_id} / PDB {pdb_id}")
            print(f"{'='*60}")

        try:
            entry_df = _process_single_entry(
                bmrb_id=bmrb_id,
                pdb_id=pdb_id,
                verbose=verbose,
                data_dir=output_dir,
            )
            if entry_df is None or entry_df.empty:
                if verbose:
                    print(f"  [{bmrb_id}/{pdb_id}] No usable data — skipping")
                results_log.append({"bmrb": bmrb_id, "pdb": pdb_id, "status": "no_data", "shifts": 0})
                continue

            # Cache the raw labeled shifts
            entry_df.to_csv(cache_file, index=False)
            all_raw_shifts.append(entry_df)
            results_log.append({
                "bmrb": bmrb_id, "pdb": pdb_id, "status": "ok",
                "shifts": len(entry_df),
                "helix": (entry_df["ss_class"] == "helix").sum(),
                "strand": (entry_df["ss_class"] == "strand").sum(),
                "coil": (entry_df["ss_class"] == "coil").sum(),
            })
            if verbose:
                print(f"  [{bmrb_id}/{pdb_id}] ✓ {len(entry_df)} shifts — "
                      f"helix={results_log[-1]['helix']} "
                      f"strand={results_log[-1]['strand']} "
                      f"coil={results_log[-1]['coil']}")

            # Be a good citizen — don't hammer BMRB/RCSB
            time.sleep(0.5)

        except Exception as e:
            if verbose:
                print(f"  [{bmrb_id}/{pdb_id}] ERROR: {e}")
            results_log.append({"bmrb": bmrb_id, "pdb": pdb_id, "status": f"error: {e}", "shifts": 0})
            continue

    # ---------- aggregate ----------
    if not all_raw_shifts:
        print("\n[BATCH] No entries processed successfully.")
        return pd.DataFrame()

    combined = pd.concat(all_raw_shifts, ignore_index=True)
    if verbose:
        print(f"\n{'='*60}")
        print(f"  BATCH AGGREGATION")
        print(f"{'='*60}")
        print(f"  Total raw shifts: {len(combined)}")
        print(f"  SS breakdown: {combined['ss_class'].value_counts().to_dict()}")
        print(f"  Unique residues: {combined['residue_name'].nunique()}")
        print(f"  Unique atoms: {combined['atom'].nunique()}")

    stats_df = _compute_merged_stats(combined)

    # Save outputs
    out_path = output_dir / "reference_db.csv"
    stats_df.to_csv(out_path, index=False)
    if verbose:
        print(f"\n  Reference DB saved → {out_path}")
        print(f"  Rows: {len(stats_df)}")

    # Save batch log
    log_df = pd.DataFrame(results_log)
    log_path = output_dir / "batch_log.csv"
    log_df.to_csv(log_path, index=False)
    if verbose:
        print(f"  Batch log saved → {log_path}")
        _print_batch_summary(log_df, stats_df)

    return stats_df


def _process_single_entry(bmrb_id, pdb_id, verbose, data_dir):
    """
    For one BMRB entry: fetch shifts, get SS labels, merge, return labeled DataFrame.

    Returns DataFrame with columns:
        bmrb_id, pdb_id, seq_id, residue_name, residue_one, atom, shift, ss_class
    """
    from bmrb_module import fetch_bmrb_star, parse_shifts
    from dssp_module import (
        download_cif, extract_ss_from_cif,
        download_pdb, align_sequences_and_map
    )
    from bmrb_module import extract_sequence_from_bmrb
    from stats_module import RESIDUE_ONE_LETTER

    # 1. Fetch BMRB
    if verbose:
        print(f"  Fetching BMRB {bmrb_id}...")
    nmrstar_text = fetch_bmrb_star(bmrb_id)
    df_shifts = parse_shifts(nmrstar_text)

    if df_shifts.empty:
        return None

    if verbose:
        print(f"  Parsed {len(df_shifts)} shifts")

    # Filter to backbone carbons + nitrogen (expand as needed)
    KEEP_ATOMS = {'CA', 'CB', 'C', 'N', 'HA', 'H'}
    df_shifts = df_shifts[df_shifts['atom'].isin(KEEP_ATOMS)].copy()
    if df_shifts.empty:
        return None

    # 2. Get SS labels from PDB/CIF
    if verbose:
        print(f"  Fetching structure {pdb_id}...")

    ss_map = {}
    try:
        cif_path = download_cif(pdb_id, data_dir=str(data_dir))
        ss_map = extract_ss_from_cif(cif_path)
        if verbose:
            counts = pd.Series(list(ss_map.values())).value_counts().to_dict()
            print(f"  SS from CIF: {counts}")
    except Exception as e:
        if verbose:
            print(f"  CIF SS failed ({e}), trying PDB HELIX/SHEET...")
        try:
            pdb_path = download_pdb(pdb_id, data_dir=str(data_dir))
            ss_map = _extract_ss_from_pdb_records(pdb_path)
        except Exception as e2:
            if verbose:
                print(f"  PDB SS also failed ({e2}). Using coil for all residues.")
            ss_map = {}

    # 3. Align BMRB → PDB residue numbering
    bmrb_seq = extract_sequence_from_bmrb(nmrstar_text)
    if bmrb_seq and ss_map:
        try:
            pdb_residues = sorted(ss_map.keys())
            seq_mapping = align_sequences_and_map(bmrb_seq, pdb_residues)
        except Exception as e:
            if verbose:
                print(f"  Alignment failed ({e}), using direct numbering")
            seq_mapping = None
    else:
        seq_mapping = None

    # 4. Assign SS to each shift row
    def get_ss(seq_id):
        if seq_mapping is not None:
            pdb_res = seq_mapping.get(int(seq_id))
            if pdb_res is not None:
                return ss_map.get(pdb_res, 'coil')
        # Direct lookup fallback
        try:
            return ss_map.get(int(seq_id), 'coil' if ss_map else 'unknown')
        except (ValueError, TypeError):
            return 'unknown'

    df_shifts['ss_class'] = df_shifts['seq_id'].apply(get_ss)
    df_shifts['bmrb_id'] = bmrb_id
    df_shifts['pdb_id'] = pdb_id

    # Add one-letter residue code
    df_shifts['residue_one'] = df_shifts['residue_name'].apply(
        lambda r: RESIDUE_ONE_LETTER.get(r.upper(), 'X') if pd.notna(r) else 'X'
    )

    return df_shifts


def _compute_merged_stats(combined_df):
    """
    Compute shift statistics aggregated across all entries.
    Groups by (residue_one, atom, ss_class) and computes descriptive stats.
    Includes entry_count so you know how many proteins contributed.
    """
    ATOM_TYPES = {'CA', 'CB', 'C', 'N', 'HA', 'H'}
    df = combined_df[combined_df['atom'].isin(ATOM_TYPES)].copy()
    df = df[df['ss_class'] != 'unknown'].copy()

    rows = []
    for (res1, atom, ss_class), group in df.groupby(['residue_one', 'atom', 'ss_class']):
        shifts = group['shift'].dropna().values
        if len(shifts) == 0:
            continue
        entry_count = group['bmrb_id'].nunique() if 'bmrb_id' in group.columns else 1
        rows.append({
            'residue': res1,
            'atom': atom,
            'ss_class': ss_class,
            'count': len(shifts),
            'entry_count': entry_count,
            'mean': float(np.mean(shifts)),
            'median': float(np.median(shifts)),
            'std': float(np.std(shifts)),
            'min': float(np.min(shifts)),
            'max': float(np.max(shifts)),
            'q25': float(np.percentile(shifts, 25)),
            'q75': float(np.percentile(shifts, 75)),
        })

    stats_df = pd.DataFrame(rows)
    if not stats_df.empty:
        stats_df = stats_df.sort_values(['residue', 'atom', 'ss_class']).reset_index(drop=True)
    return stats_df


def _extract_ss_from_pdb_records(pdb_path):
    """
    Fallback: read HELIX/SHEET records from legacy PDB file.
    Returns {residue_number: 'helix'|'strand'|'coil'}.
    """
    helix_ranges = []
    sheet_ranges = []

    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith('HELIX'):
                try:
                    start = int(line[21:25].strip())
                    end   = int(line[33:37].strip())
                    helix_ranges.append((start, end))
                except ValueError:
                    pass
            elif line.startswith('SHEET'):
                try:
                    start = int(line[22:26].strip())
                    end   = int(line[33:37].strip())
                    sheet_ranges.append((start, end))
                except ValueError:
                    pass

    # Get all residue numbers from ATOM records
    residue_nums = set()
    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                try:
                    residue_nums.add(int(line[22:26].strip()))
                except ValueError:
                    pass

    ss_map = {}
    for res in residue_nums:
        if any(s <= res <= e for s, e in helix_ranges):
            ss_map[res] = 'helix'
        elif any(s <= res <= e for s, e in sheet_ranges):
            ss_map[res] = 'strand'
        else:
            ss_map[res] = 'coil'
    return ss_map


def _print_batch_summary(log_df, stats_df):
    """Print a clean batch run summary."""
    ok = log_df[log_df['status'] == 'ok']
    cached = log_df[log_df['status'] == 'cached']
    skipped = log_df[log_df['status'] == 'skipped']
    errors = log_df[~log_df['status'].isin(['ok', 'cached', 'skipped'])]

    total_shifts = log_df[log_df['status'].isin(['ok', 'cached'])]['shifts'].sum()

    print(f"\n{'='*60}")
    print(f"  BATCH SUMMARY")
    print(f"{'='*60}")
    print(f"  Processed:  {len(ok)} entries (fresh)")
    print(f"  Cached:     {len(cached)} entries")
    print(f"  Skipped:    {len(skipped)} entries")
    print(f"  Errors:     {len(errors)} entries")
    print(f"  Total shifts in reference DB: {total_shifts}")

    if not stats_df.empty and 'ss_class' in stats_df.columns:
        ss_counts = stats_df.groupby('ss_class')['count'].sum()
        print(f"\n  Reference DB SS breakdown:")
        for ss, cnt in ss_counts.items():
            print(f"    {ss:8s}: {cnt} observations")

        # Show secondary chemical shift effect if we have helix+strand+coil
        print(f"\n  Secondary chemical shift validation (Ala CA, if present):")
        ala_ca = stats_df[(stats_df['residue'] == 'A') & (stats_df['atom'] == 'CA')]
        for _, row in ala_ca.iterrows():
            print(f"    Ala CA [{row['ss_class']:6s}]: mean={row['mean']:.2f} ppm  n={row['count']}")

    if len(errors) > 0:
        print(f"\n  Failed entries:")
        for _, row in errors.iterrows():
            print(f"    BMRB {row['bmrb']} / PDB {row['pdb']}: {row['status']}")

    print(f"{'='*60}")


# ---------------------------------------------------------------------------
# Convenience function: load existing reference DB
# ---------------------------------------------------------------------------

def load_reference_db(data_dir=None):
    """
    Load the merged reference database from disk.
    Returns None if not found (run run_batch() first).
    """
    if data_dir is None:
        data_dir = Path(__file__).resolve().parent.parent / "data"
    path = Path(data_dir) / "reference_db.csv"
    if not path.exists():
        print(f"[BATCH] Reference DB not found at {path}")
        print("[BATCH] Run: python src/pipeline.py --batch")
        return None
    df = pd.read_csv(path)
    print(f"[BATCH] Loaded reference DB: {len(df)} rows from {path}")
    return df


def get_entry_list():
    """Return the curated entry list as a DataFrame for inspection."""
    rows = []
    for entry in SOLID_STATE_ENTRIES:
        bmrb_id, pdb_id, desc = entry[0], entry[1], entry[2]
        in_skip = bmrb_id in SKIP_ENTRIES
        rows.append({
            "bmrb_id": bmrb_id,
            "pdb_id": pdb_id,
            "description": desc,
            "skip": in_skip,
        })
    return pd.DataFrame(rows)
