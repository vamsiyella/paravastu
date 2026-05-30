"""
pipeline.py — Main orchestrator for the Paravastu NMR structural annotation pipeline.

Usage:
    python src/pipeline.py --bmrb 17561 --pdb 2LBH          # single entry (Phase 1/2)
    python src/pipeline.py --batch                            # Phase 3: batch all curated entries
    python src/pipeline.py --batch --entries 17561:2LBH 15409:2KIB  # batch with custom pairs
    python src/pipeline.py --list-entries                     # show curated entry list
    python src/pipeline.py --csv data/my_shifts.csv --pdb 2LBH      # from CSV
    python src/pipeline.py --bmrb 17561 --no-dssp            # skip DSSP
"""

import sys
import argparse
import pandas as pd
from pathlib import Path

# ── make src/ importable regardless of working directory ──────────────────
SRC_DIR = Path(__file__).resolve().parent
ROOT_DIR = SRC_DIR.parent
DATA_DIR = ROOT_DIR / "data"
RESULTS_DIR = ROOT_DIR / "results"
sys.path.insert(0, str(SRC_DIR))

from bmrb_module import (
    fetch_bmrb_text,
    load_bmrb_from_file,
    parse_shifts,
    extract_sequence,
    extract_pdb_id,
    compute_coverage,
    compute_shift_stats,
)
from dssp_module import (
    download_pdb,
    extract_dssp,
    extract_dssp_full,
    extract_ss_from_pdb_records,
    build_ss_segments,
    build_prediction_ss_map,
    align_and_map,
    get_pdb_sequence_from_file,
    find_mkdssp,
)
from stats_module import (
    compute_shift_stats as stats_compute,
    add_random_coil_deviation,
    ShiftPredictor,
)
from viz_module import generate_all_plots
from viz_module_phase4 import generate_phase4_plots
from ml_module import run_ml_pipeline
from assignment_engine import run_assignment_engine


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def print_section(title: str):
    print(f"\n{'='*60}")
    print(f"  {title}")
    print('='*60)


def merge_shifts_with_dssp(
    shifts_df: pd.DataFrame,
    dssp_map: dict,
    seq_mapping: dict = None,
) -> pd.DataFrame:
    """
    Join chemical shifts with DSSP secondary structure labels.
    seq_mapping: {bmrb_seq_id -> pdb_residue_number} (from alignment)
    If None, assumes 1:1 numbering.

    If dssp_map is empty AND shifts_df already has a non-unknown ss_class column,
    those existing labels are preserved rather than overwritten.
    """
    df = shifts_df.copy()

    if not dssp_map:
        if 'ss_class' in df.columns and df['ss_class'].nunique() > 1:
            print("[MERGE] No DSSP map supplied — keeping existing ss_class labels from CSV.")
            return df
        df['ss_class'] = 'unknown'
        return df

    def get_ss(seq_id):
        if seq_mapping is not None:
            pdb_res = seq_mapping.get(seq_id)
            if pdb_res is None:
                return 'unknown'
            return dssp_map.get(pdb_res, 'unknown')
        return dssp_map.get(seq_id, 'unknown')

    df['ss_class'] = df['seq_id'].apply(get_ss)
    return df


# ---------------------------------------------------------------------------
# Main pipeline (single entry)
# ---------------------------------------------------------------------------

def run_pipeline(
    bmrb_id: int = 17561,
    pdb_id: str = "2LBH",
    run_dssp: bool = True,
    save_results: bool = True,
) -> dict:
    """
    Full pipeline run for a single BMRB entry. Returns dict with all results.
    """
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    results = {}

    # ── Step 1: Load BMRB data ─────────────────────────────────────────────
    print_section(f"STEP 1 — Loading BMRB {bmrb_id}")

    local_str = DATA_DIR / f"bmr{bmrb_id}_3.str"
    nmrstar_text = None

    if local_str.exists():
        raw = local_str.read_text(errors="replace")
        if raw[:20].strip().startswith("data_"):
            print(f"Using local BMRB file: {local_str}")
            nmrstar_text = raw
        else:
            print(f"Cached file invalid (starts: {raw[:30]!r}). Deleting and re-fetching.")
            local_str.unlink()

    if nmrstar_text is None:
        try:
            nmrstar_text = fetch_bmrb_text(bmrb_id, cache_dir=DATA_DIR)
        except RuntimeError as e:
            print(f"ERROR: {e}")
            print(f"  Manual download: https://bmrb.io/ftp/pub/bmrb/entry_directories/bmr{bmrb_id}/bmr{bmrb_id}_3.str")
            print(f"  Save to: {local_str}")
            return {}

    # ── Step 2: Parse shifts + sequence ───────────────────────────────────
    print_section("STEP 2 — Parsing shifts and sequence")

    shifts_df = parse_shifts(nmrstar_text)
    print(f"Parsed {len(shifts_df)} shift records")
    print(shifts_df.head(5).to_string())

    sequence = extract_sequence(nmrstar_text)
    if sequence:
        print(f"\nSequence ({len(sequence)} aa):\n{sequence}")
    else:
        print("WARNING: Could not extract sequence")

    if pdb_id is None:
        pdb_id = extract_pdb_id(nmrstar_text)
        print(f"PDB ID from BMRB: {pdb_id}")

    results['shifts_df'] = shifts_df
    results['sequence']  = sequence
    results['pdb_id']    = pdb_id

    # ── Step 3: Coverage analysis ──────────────────────────────────────────
    print_section("STEP 3 — Coverage analysis")

    coverage_df = compute_coverage(shifts_df)
    print(f"Residue coverage (first 10):")
    print(coverage_df.head(10).to_string())
    print(f"\nAverage atoms observed per residue: {coverage_df['n_atoms'].mean():.2f}")
    results['coverage_df'] = coverage_df

    if save_results:
        coverage_df.to_csv(RESULTS_DIR / f"coverage_{bmrb_id}.csv", index=False)

    # ── Step 4: DSSP ──────────────────────────────────────────────────────
    dssp_map = {}
    seq_mapping = None

    if run_dssp:
        print_section(f"STEP 4 — DSSP from PDB {pdb_id}")
        mkdssp = find_mkdssp()
        if mkdssp is None:
            print("WARNING: mkdssp not found. Skipping DSSP.")
            print("Fix: conda install -c salilab dssp")
            run_dssp = False

    if run_dssp and pdb_id:
        pdb_path = download_pdb(pdb_id, DATA_DIR)

        dssp_df = None
        try:
            dssp_df = extract_dssp_full(pdb_path)
            print(f"[DSSP] mkdssp succeeded.")
        except Exception as e:
            print(f"[DSSP] mkdssp failed: {e}")
            print("[DSSP] Falling back to CIF HELIX/SHEET records...")
            try:
                dssp_df = extract_ss_from_pdb_records(pdb_path)
                print("[DSSP] CIF fallback succeeded.")
            except Exception as e2:
                print(f"[DSSP] Fallback also failed: {e2}")
                print("Continuing without structural labels.")

        if dssp_df is not None and not dssp_df.empty:
            print(f"\nSS assignments (first 10):")
            print(dssp_df.head(10).to_string())
            print(f"\nSS distribution:\n{dssp_df['ss_class'].value_counts()}")

            dssp_map = dict(zip(dssp_df['residue_number'], dssp_df['ss_class']))
            results['dssp_df'] = dssp_df

            segments = build_ss_segments(dssp_map)
            print(f"\nSS segments:")
            for seg in segments:
                print(f"  {seg['ss_class']:7s} res {seg['start']:3d}-{seg['end']:3d}  ({seg['length']} aa)")
            results['segments'] = segments

            if sequence:
                pdb_seq = get_pdb_sequence_from_file(pdb_path)
                print(f"\nBMRB seq ({len(sequence)} aa): {sequence[:40]}...")
                print(f"PDB  seq ({len(pdb_seq)} aa):  {pdb_seq[:40]}...")
                seq_mapping = align_and_map(sequence, pdb_seq)
                results['seq_mapping'] = seq_mapping

            if save_results:
                dssp_df.to_csv(RESULTS_DIR / f"dssp_{pdb_id}.csv", index=False)

    # ── Step 5: Merge shifts + DSSP ───────────────────────────────────────
    print_section("STEP 5 — Merging shifts with structural labels")

    merged_df = merge_shifts_with_dssp(shifts_df, dssp_map, seq_mapping)
    ss_counts = merged_df['ss_class'].value_counts()
    print(f"Shift records by SS class:\n{ss_counts}")
    results['merged_df'] = merged_df

    if save_results:
        merged_df.to_csv(RESULTS_DIR / f"merged_shifts_{bmrb_id}.csv", index=False)
        print(f"Saved merged shifts → results/merged_shifts_{bmrb_id}.csv")

    # ── Step 6: Compute statistics ─────────────────────────────────────────
    print_section("STEP 6 — Computing shift statistics by SS class")

    stats_df = stats_compute(merged_df)
    stats_df = add_random_coil_deviation(stats_df)
    print(f"Statistics table ({len(stats_df)} rows):")
    print(stats_df.head(20).to_string())
    results['stats_df'] = stats_df

    if save_results:
        stats_df.to_csv(RESULTS_DIR / f"stats_{bmrb_id}.csv", index=False)
        print(f"Saved statistics → results/stats_{bmrb_id}.csv")

    # ── Step 7: Prediction demo ────────────────────────────────────────────
    if sequence and stats_df is not None and not stats_df.empty:
        print_section("STEP 7 — Shift prediction demo")

        predictor = ShiftPredictor(stats_df)
        seq_mapping_local = results.get('seq_mapping')
        dssp_df_local     = results.get('dssp_df')
        if dssp_map and dssp_df_local is not None and seq_mapping_local:
            ss_map_for_pred = build_prediction_ss_map(
                dssp_df_local, seq_mapping_local, len(sequence)
            )
        elif dssp_map:
            ss_map_for_pred = {i+1: dssp_map.get(i+1, 'coil') for i in range(len(sequence))}
        else:
            ss_map_for_pred = {i+1: 'coil' for i in range(len(sequence))}

        print(f"Prediction SS map: { {k: list(ss_map_for_pred.values()).count(k) for k in ['coil','strand','helix'] if k in ss_map_for_pred.values()} }")
        predictions = predictor.predict(sequence, ss_map_for_pred)
        predictions = add_random_coil_deviation(predictions.rename(columns={'predicted_shift': 'shift'}))
        predictions.rename(columns={'shift': 'predicted_shift'}, inplace=True)

        print(f"Sample predictions:")
        print(predictions.head(15).to_string())
        results['predictions'] = predictions

        if save_results:
            predictions.to_csv(RESULTS_DIR / f"predictions_{bmrb_id}.csv", index=False)
            print(f"Saved predictions → results/predictions_{bmrb_id}.csv")

    # ── Step 8: Visualizations ────────────────────────────────────────────
    print_section("STEP 8 — Generating plots")
    try:
        results['bmrb_id'] = bmrb_id
        generate_all_plots(results, output_dir=RESULTS_DIR)
    except Exception as e:
        print(f"[VIZ] Plot generation failed: {e}")

    print_section("STEP 8b — Phase 4 spectrum simulation")
    try:
        generate_phase4_plots(results, output_dir=RESULTS_DIR, bmrb_id=str(bmrb_id))
    except Exception as e:
        print(f"[Phase 4] {e}")

    # ── Step 9: ML model training ──────────────────────────────────────────
    if merged_df is not None and not merged_df.empty:
        labeled = merged_df[merged_df['ss_class'].isin(['helix', 'strand', 'coil'])]
        if len(labeled) >= 15 and labeled['ss_class'].nunique() >= 2:
            print_section("STEP 9 — Training ML models")
            try:
                ml_results = run_ml_pipeline(merged_df, results_dir=RESULTS_DIR)
                results['ml'] = ml_results
                if ml_results:
                    print(f"  RF  CV: {ml_results['rf']['cv_scores'].mean():.3f}")
                    print(f"  XGB CV: {ml_results['xgb']['cv_scores'].mean():.3f}")
                    print(f"  Top features: {', '.join(ml_results['rf']['feature_importance'].head(5).index.tolist())}")
            except Exception as e:
                print(f"[ML] Training failed: {e}")
        else:
            print_section("STEP 9 — ML skipped (not enough labeled samples)")
            print(f"  Run --batch to build a larger training set.")

    # ── Step 10: Phase 5 Assignment Engine ───────────────────────────────────
    print_section("STEP 10 — Assignment Engine (Phase 5)")
    if stats_df is not None and sequence and dssp_map:
        # Placeholder: replace with real experimental peaks when you have them
        # Format: {'CA': [list of ppm values], 'CB': [...], 'N': [...]}
        # ca_cb_pairs: [(ca_ppm, cb_ppm), ...] for 2D joint assignment
        print("  [Phase 5] Pass experimental_peaks dict to run_assignment_engine()")
        print("  Example:")
        print("    from assignment_engine import run_assignment_engine")
        print("    run_assignment_engine({'CA': [55.3, 62.1, 45.0, ...]}, stats_df, sequence, ss_map, RESULTS_DIR)")

        # ── Summary ───────────────────────────────────────────────────────────
        print_section("PIPELINE COMPLETE")
        print(f"  BMRB entry:      {bmrb_id}")
        print(f"  PDB structure:   {pdb_id or 'none'}")
        print(f"  Shift records:   {len(shifts_df)}")
        print(f"  Sequence length: {len(sequence) if sequence else 'unknown'}")
        print(f"  DSSP residues:   {len(dssp_map)}")
        print(f"  Stats rows:      {len(stats_df) if stats_df is not None else 0}")
        print(f"  Plots:    {len(list(RESULTS_DIR.glob('*.png')))} PNG files")
        print(f"  CSVs:     {len(list(RESULTS_DIR.glob('*.csv')))} CSV files")
        print(f"  Models:   {len(list(RESULTS_DIR.glob('*.joblib')))} trained models")
        print(f"\nOutputs saved to: {RESULTS_DIR}")

        return results


# ---------------------------------------------------------------------------
# Phase 3: Batch processing
# ---------------------------------------------------------------------------

# Curated solid-state NMR BMRB entries with linked PDB structures.
# Selection criteria: ssNMR, 13C/15N shifts deposited, PDB structure available.
SOLID_STATE_ENTRIES = [
    # (bmrb_id, pdb_id, description)
    (17561, "2LBH", "EETI-II knottin — beta-sheet + coil (Phase 1/2 entry)"),
    (15409, "2KIB", "HET-s prion domain — beta-solenoid"),
    (16318, "2JWU", "ubiquitin microcrystals — alpha/beta mixed"),
    (15380, "1LY2", "GB1 protein — alpha/beta mixed"),
    (15865, "1M8M", "SH3 domain — all beta"),
    (6838,  "1YMZ", "fd coat protein — helix-rich"),
    (17557, "2KSF", "M2 proton channel — helix bundle"),
    (16299, "2JSV", "thioredoxin microcrystals — alpha/beta")
    #(5969,  "1H4W", "BPTI — disulfide-rich beta"),
    #(17948, "2NUZ", "calmodulin — helix-rich"),
]


def run_pipeline_batch(
    custom_pairs=None,
    no_ml: bool = False,
    verbose: bool = True,
) -> pd.DataFrame:
    """
    Phase 3: Process multiple BMRB entries and build a merged reference database.

    Parameters
    ----------
    custom_pairs : list of (bmrb_id, pdb_id) or None
        If None, uses the built-in SOLID_STATE_ENTRIES list.
        Pass your own: [(17561, '2LBH'), (15409, '2KIB')]
    no_ml : bool
        Skip ML re-training after batch (faster when just rebuilding the DB).
    verbose : bool
        Print per-entry progress.

    Returns
    -------
    pd.DataFrame  — merged reference_db.csv
    """
    import time

    DATA_DIR.mkdir(parents=True, exist_ok=True)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    cache_dir = DATA_DIR / "batch_cache"
    cache_dir.mkdir(exist_ok=True)

    pairs = custom_pairs or [(e[0], e[1]) for e in SOLID_STATE_ENTRIES]

    print_section("PHASE 3 — BATCH REFERENCE DATABASE BUILD")
    print(f"  Entries to process: {len(pairs)}")
    print(f"  Cache directory:    {cache_dir}")
    print(f"  Output:             {DATA_DIR / 'reference_db.csv'}")

    all_raw = []
    log = []

    for bmrb_id, pdb_id in pairs:
        cache_file = cache_dir / f"bmr{bmrb_id}_{pdb_id}_raw.csv"

        # Use cached raw shifts if available
        if cache_file.exists():
            if verbose:
                print(f"\n  [{bmrb_id}/{pdb_id}] Loading from cache...")
            try:
                df = pd.read_csv(cache_file)
                all_raw.append(df)
                ss = df['ss_class'].value_counts().to_dict()
                log.append({"bmrb": bmrb_id, "pdb": pdb_id, "status": "cached",
                             "shifts": len(df), **ss})
                if verbose:
                    print(f"    {len(df)} shifts — {ss}")
                continue
            except Exception as e:
                if verbose:
                    print(f"    Cache read failed ({e}), re-processing...")

        if verbose:
            print(f"\n  [{bmrb_id}/{pdb_id}] Processing...")

        try:
            res = run_pipeline(
                bmrb_id=bmrb_id,
                pdb_id=pdb_id,
                run_dssp=True,
                save_results=False,
            )

            merged = res.get('merged_df')
            if merged is None or merged.empty:
                if verbose:
                    print(f"    No usable data.")
                log.append({"bmrb": bmrb_id, "pdb": pdb_id, "status": "no_data", "shifts": 0})
                continue

            # Tag with source identifiers and cache
            merged['bmrb_id'] = bmrb_id
            merged['pdb_id']  = pdb_id
            merged.to_csv(cache_file, index=False)
            all_raw.append(merged)

            ss = merged['ss_class'].value_counts().to_dict()
            log.append({"bmrb": bmrb_id, "pdb": pdb_id, "status": "ok",
                         "shifts": len(merged), **ss})
            if verbose:
                print(f"    ✓ {len(merged)} shifts — {ss}")

            time.sleep(0.3)  # be polite to BMRB/RCSB

        except Exception as e:
            if verbose:
                print(f"    ERROR: {e}")
            log.append({"bmrb": bmrb_id, "pdb": pdb_id, "status": f"error: {e}", "shifts": 0})
            continue

    # ── Aggregate ────────────────────────────────────────────────────────
    if not all_raw:
        print("\n[BATCH] No entries processed successfully.")
        return pd.DataFrame()

    combined = pd.concat(all_raw, ignore_index=True)

    print_section("BATCH AGGREGATION")
    print(f"  Total raw shifts: {len(combined)}")
    print(f"  SS breakdown:     {combined['ss_class'].value_counts().to_dict()}")

    # Build merged reference stats
    reference_df = _compute_reference_stats(combined)
    ref_path = DATA_DIR / "reference_db.csv"
    reference_df.to_csv(ref_path, index=False)
    print(f"  Reference DB:     {len(reference_df)} rows → {ref_path}")

    # Save batch log
    log_df = pd.DataFrame(log)
    log_df.to_csv(DATA_DIR / "batch_log.csv", index=False)

    # Print summary
    ok = log_df[log_df['status'].isin(['ok', 'cached'])]
    err = log_df[~log_df['status'].isin(['ok', 'cached', 'no_data'])]
    print(f"\n  Processed: {len(ok)} entries")
    if len(err):
        print(f"  Errors:    {len(err)} entries")
        for _, row in err.iterrows():
            print(f"    BMRB {row['bmrb']} / {row['pdb']}: {row['status']}")

    # Validate secondary chemical shift effect
    ala_ca = reference_df[(reference_df['residue'] == 'A') & (reference_df['atom'] == 'CA')]
    if len(ala_ca) > 1:
        print(f"\n  Secondary chemical shift check (Ala CA):")
        for _, row in ala_ca.iterrows():
            print(f"    [{row['ss_class']:6s}] mean={row['mean']:.2f} ppm  n={row['count']}")

    # ── ML training on full combined dataset ─────────────────────────────
    if not no_ml:
        labeled = combined[combined['ss_class'].isin(['helix', 'strand', 'coil'])]
        if len(labeled) >= 30 and labeled['ss_class'].nunique() >= 2:
            print_section("PHASE 3 — ML TRAINING ON REFERENCE DB")
            print(f"  Training on {len(labeled)} labeled shifts from {len(ok)} entries")
            try:
                ml_results = run_ml_pipeline(combined, results_dir=RESULTS_DIR)
                if ml_results:
                    print(f"\n  RF  CV accuracy: {ml_results['rf']['cv_scores'].mean():.1%}")
                    print(f"  XGB CV accuracy: {ml_results['xgb']['cv_scores'].mean():.1%}")
                    print(f"  Models saved → {RESULTS_DIR}/")
            except Exception as e:
                print(f"  ML training failed: {e}")
        else:
            print("\n[BATCH] Not enough labeled data for ML training.")
            print("  Need at least 30 labeled shifts across 2+ SS classes.")

    print_section("PHASE 3 COMPLETE")
    print(f"  Reference DB:  {ref_path}")
    print(f"  Batch log:     {DATA_DIR / 'batch_log.csv'}")
    print(f"  Cached shifts: {cache_dir}/")
    print(f"\n  Next: run --bmrb 17561 --pdb 2LBH to use the new reference DB for predictions.")

    return reference_df


def _compute_reference_stats(combined_df: pd.DataFrame) -> pd.DataFrame:
    """
    Aggregate shift statistics across all entries in the batch.
    Groups by (residue, atom, ss_class). Adds entry_count column.
    Only includes backbone atoms most useful for structure analysis.
    """
    import numpy as np

    KEEP_ATOMS = {'CA', 'CB', 'C', 'N', 'H', 'HA'}
    df = combined_df[combined_df['atom'].isin(KEEP_ATOMS)].copy()
    df = df[df['ss_class'].isin(['helix', 'strand', 'coil'])].copy()

    rows = []
    for (res, atom, ss), group in df.groupby(['residue', 'atom', 'ss_class']):
        shifts = group['shift'].dropna().values
        if len(shifts) == 0:
            continue
        entry_count = group['bmrb_id'].nunique() if 'bmrb_id' in group.columns else 1
        rows.append({
            'residue':     res,
            'atom':        atom,
            'ss_class':    ss,
            'count':       len(shifts),
            'entry_count': int(entry_count),
            'mean':        float(np.mean(shifts)),
            'median':      float(np.median(shifts)),
            'std':         float(np.std(shifts)),
            'min':         float(np.min(shifts)),
            'max':         float(np.max(shifts)),
        })

    stats_df = pd.DataFrame(rows)
    if not stats_df.empty:
        stats_df = add_random_coil_deviation(stats_df.rename(columns={'mean': 'shift'}))
        stats_df = stats_df.rename(columns={'shift': 'mean'})
        stats_df = stats_df.sort_values(['residue', 'atom', 'ss_class']).reset_index(drop=True)
    return stats_df


# ---------------------------------------------------------------------------
# CSV-based entry point
# ---------------------------------------------------------------------------

def run_pipeline_from_csv(
    csv_path: str,
    pdb_id: str = None,
    sequence: str = None,
    run_dssp: bool = True,
    save_results: bool = True,
    bmrb_id: int = 0,
) -> dict:
    """
    Run the pipeline starting from an already-parsed shift CSV.

    Accepts two CSV formats:
    - STATS CSV:  residue, atom, ss_class, count, mean, median, std, min, max
    - RAW CSV:    seq_id, residue, atom, shift, ss_class
    """
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    results = {'bmrb_id': bmrb_id}

    print_section(f"Loading shifts from CSV: {csv_path}")
    shifts_df = pd.read_csv(csv_path)

    COARSE = {'H': 'helix', 'G': 'helix', 'I': 'helix',
              'E': 'strand', 'B': 'strand',
              'C': 'coil', 'T': 'coil', 'S': 'coil', '-': 'coil', ' ': 'coil',
              'helix': 'helix', 'strand': 'strand', 'coil': 'coil',
              'unknown': 'unknown'}

    if 'ss_class' in shifts_df.columns:
        shifts_df['ss_class'] = shifts_df['ss_class'].map(COARSE).fillna('unknown')
    else:
        shifts_df['ss_class'] = 'unknown'

    if 'residue' not in shifts_df.columns and 'residue_name' in shifts_df.columns:
        THREE_TO_ONE = {
            'ALA':'A','CYS':'C','ASP':'D','GLU':'E','PHE':'F','GLY':'G',
            'HIS':'H','ILE':'I','LYS':'K','LEU':'L','MET':'M','ASN':'N',
            'PRO':'P','GLN':'Q','ARG':'R','SER':'S','THR':'T','VAL':'V',
            'TRP':'W','TYR':'Y',
        }
        shifts_df['residue'] = shifts_df['residue_name'].apply(
            lambda x: THREE_TO_ONE.get(str(x).upper(), 'X')
        )

    lbl = bmrb_id or Path(csv_path).stem
    is_stats_csv = 'mean' in shifts_df.columns and 'seq_id' not in shifts_df.columns

    if is_stats_csv:
        print("Detected: STATS CSV (pre-aggregated)")
        print(f"  Rows: {len(shifts_df)}, atoms: {sorted(shifts_df['atom'].unique())}")
        stats_df = shifts_df.copy()
        stats_df['ss_class'] = stats_df['ss_class'].map(COARSE).fillna('unknown')
        merged_df = None
        coverage_df = None
    else:
        print("Detected: RAW SHIFTS CSV (one row per observed shift)")
        shifts_df['ss_class'] = shifts_df['ss_class'].map(COARSE).fillna('unknown')
        print(f"  Rows: {len(shifts_df)}, SS: {shifts_df['ss_class'].unique().tolist()}")
        stats_df = None
        merged_df = None
        coverage_df = None

    results['shifts_df'] = shifts_df
    results['sequence']  = sequence
    results['pdb_id']    = pdb_id

    if not is_stats_csv:
        print_section("Coverage analysis")
        coverage_df = compute_coverage(shifts_df)
        print(f"Residues: {len(coverage_df)}, avg atoms/res: {coverage_df['n_atoms'].mean():.2f}")
        results['coverage_df'] = coverage_df
        if save_results:
            coverage_df.to_csv(RESULTS_DIR / f"coverage_{lbl}.csv", index=False)

    # ── DSSP ──────────────────────────────────────────────────────────────
    dssp_map = {}
    seq_mapping = None

    if run_dssp and pdb_id:
        print_section(f"DSSP from PDB {pdb_id}")
        mkdssp = find_mkdssp()
        if mkdssp is None:
            print("WARNING: mkdssp not found. Run: conda install -c salilab dssp")
        else:
            pdb_path_dl = download_pdb(pdb_id, DATA_DIR)
            dssp_df = None
            try:
                dssp_df = extract_dssp_full(pdb_path_dl)
                print(f"[DSSP] mkdssp succeeded.")
            except Exception as e:
                print(f"[DSSP] mkdssp failed: {e}")
                print("[DSSP] Falling back to CIF records...")
                try:
                    dssp_df = extract_ss_from_pdb_records(pdb_path_dl)
                    print("[DSSP] CIF fallback succeeded.")
                except Exception as e2:
                    print(f"[DSSP] Fallback also failed: {e2}")

            if dssp_df is not None and not dssp_df.empty:
                dssp_map = dict(zip(dssp_df['residue_number'], dssp_df['ss_class']))
                segments = build_ss_segments(dssp_map)
                results['dssp_df']  = dssp_df
                results['segments'] = segments
                print(f"SS labels: {len(dssp_map)} residues")
                print(dssp_df['ss_class'].value_counts().to_string())
                if sequence:
                    pdb_seq     = get_pdb_sequence_from_file(pdb_path_dl)
                    seq_mapping = align_and_map(sequence, pdb_seq)
                if save_results:
                    dssp_df.to_csv(RESULTS_DIR / f"dssp_{pdb_id}.csv", index=False)

    # ── Merge / stats ─────────────────────────────────────────────────────
    if not is_stats_csv:
        print_section("Merging shifts with structural labels")
        merged_df = merge_shifts_with_dssp(shifts_df, dssp_map, seq_mapping)
        results['merged_df'] = merged_df
        if save_results:
            merged_df.to_csv(RESULTS_DIR / f"merged_shifts_{lbl}.csv", index=False)
        print(merged_df['ss_class'].value_counts().to_string())

        print_section("Computing shift statistics by SS class")
        stats_df = stats_compute(merged_df)
        stats_df = add_random_coil_deviation(stats_df)
        results['stats_df'] = stats_df
        if save_results:
            stats_df.to_csv(RESULTS_DIR / f"stats_{lbl}.csv", index=False)
        print(f"{len(stats_df)} stat rows")
    else:
        print_section("Using pre-computed statistics")
        stats_df = add_random_coil_deviation(stats_df)
        results['stats_df'] = stats_df
        if save_results:
            stats_df.to_csv(RESULTS_DIR / f"stats_{lbl}.csv", index=False)
        print(f"{len(stats_df)} stat rows")
        if dssp_map:
            print(f"\nNOTE: DSSP ran. For SS-aware stats, run:")
            print(f"  python src/pipeline.py --bmrb {bmrb_id or 17561} --pdb {pdb_id}")

    # ── Prediction ────────────────────────────────────────────────────────
    if sequence and stats_df is not None and not stats_df.empty:
        print_section("Shift prediction")
        predictor = ShiftPredictor(stats_df)
        if dssp_map and results.get('dssp_df') is not None and seq_mapping:
            ss_map_for_pred = build_prediction_ss_map(
                results['dssp_df'], seq_mapping, len(sequence)
            )
        elif dssp_map:
            ss_map_for_pred = {i+1: dssp_map.get(i+1, 'coil') for i in range(len(sequence))}
        else:
            ss_map_for_pred = {i+1: 'unknown' for i in range(len(sequence))}
        preds = predictor.predict(sequence, ss_map_for_pred, atoms=['CA', 'N'])
        preds = preds.rename(columns={'predicted_shift': 'shift'})
        preds = add_random_coil_deviation(preds)
        preds.rename(columns={'shift': 'predicted_shift'}, inplace=True)
        results['predictions'] = preds
        if save_results:
            preds.to_csv(RESULTS_DIR / f"predictions_{lbl}.csv", index=False)
        print(f"Predicted shifts for {len(preds)} (residue, atom) pairs")
        print(preds.head(10).to_string())

    # ── Visualizations ────────────────────────────────────────────────────
    print_section("Generating plots")
    try:
        results['bmrb_id'] = lbl
        generate_all_plots(results, output_dir=RESULTS_DIR)
    except Exception as e:
        print(f"[VIZ] Plot generation failed: {e}")

    # ── ML ────────────────────────────────────────────────────────────────
    if not is_stats_csv and merged_df is not None:
        labeled = merged_df[merged_df['ss_class'].isin(['helix', 'strand', 'coil'])]
        if len(labeled) >= 15 and labeled['ss_class'].nunique() >= 2:
            print_section("Training ML models")
            try:
                ml_results = run_ml_pipeline(merged_df, results_dir=RESULTS_DIR)
                results['ml'] = ml_results
            except Exception as e:
                print(f"[ML] Training failed: {e}")

    print_section("PIPELINE COMPLETE")
    print(f"  Source CSV:    {csv_path}")
    print(f"  CSV type:      {'stats (pre-aggregated)' if is_stats_csv else 'raw shifts'}")
    print(f"  PDB:           {pdb_id or 'none'}")
    print(f"  DSSP residues: {len(dssp_map)}")
    print(f"  Stat rows:     {len(stats_df) if stats_df is not None else 0}")
    print(f"  Plots:   {len(list(RESULTS_DIR.glob('*.png')))} PNG files")
    print(f"  Results → {RESULTS_DIR}")

    return results


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Paravastu NMR Structural Annotation Pipeline"
    )

    # Single-entry args
    parser.add_argument("--bmrb",    type=int, default=17561, help="BMRB entry ID")
    parser.add_argument("--pdb",     type=str, default="2LBH", help="PDB ID (or 'auto')")
    parser.add_argument("--csv",     type=str, default=None,   help="Path to shifts CSV")
    parser.add_argument("--seq",     type=str, default=None,   help="Protein sequence (one-letter)")
    parser.add_argument("--no-dssp", action="store_true",      help="Skip DSSP extraction")
    parser.add_argument("--no-save", action="store_true",      help="Don't write output files")
    parser.add_argument("--phase4", action="store_true",       help="Run Phase 4 enhanced spectrum simulation (Voigt, 2D CA-CB, NMRPipe export)")
    parser.add_argument("--assign-peaks", type=str, default=None, help="Path to CSV with experimental peaks (columns: atom, shift) for Phase 5 assignment")
    # Phase 3 batch args
    batch_group = parser.add_argument_group("Phase 3 — Batch processing")
    batch_group.add_argument(
        "--batch", action="store_true",
        help="Run Phase 3: process curated BMRB entries and build reference DB"
    )
    batch_group.add_argument(
        "--entries", nargs="+", metavar="BMRB:PDB",
        help="Custom entry pairs for batch, e.g. --entries 17561:2LBH 15409:2KIB"
    )
    batch_group.add_argument(
        "--list-entries", action="store_true",
        help="List all curated solid-state BMRB entries and exit"
    )
    batch_group.add_argument(
        "--no-ml", action="store_true",
        help="Skip ML training after batch (faster)"
    )

    args = parser.parse_args()

    # ── Dispatch ────────────────────────────────────────────────────────
    if args.list_entries:
        print("\nCurated solid-state NMR entries:")
        print(f"  {'BMRB':>6}  {'PDB':<6}  Description")
        print(f"  {'-'*6}  {'-'*6}  {'-'*45}")
        for entry in SOLID_STATE_ENTRIES:
            bmrb_id, pdb_id, desc = entry
            print(f"  {bmrb_id:>6}  {pdb_id:<6}  {desc}")
        print(f"\n  Total: {len(SOLID_STATE_ENTRIES)} entries")
        sys.exit(0)

    if args.batch or args.entries:
        custom_pairs = None
        if args.entries:
            custom_pairs = []
            for pair in args.entries:
                parts = pair.split(":")
                if len(parts) != 2:
                    print(f"Error: --entries format is BMRB_ID:PDB_ID, got '{pair}'")
                    sys.exit(1)
                try:
                    custom_pairs.append((int(parts[0]), parts[1].upper()))
                except ValueError:
                    print(f"Error: BMRB ID must be an integer, got '{parts[0]}'")
                    sys.exit(1)
        run_pipeline_batch(
            custom_pairs=custom_pairs,
            no_ml=args.no_ml,
        )
        sys.exit(0)

    pdb_id = None if args.pdb == 'auto' else args.pdb

    if args.csv:
        run_pipeline_from_csv(
            csv_path=args.csv,
            pdb_id=pdb_id,
            sequence=args.seq,
            run_dssp=not args.no_dssp,
            save_results=not args.no_save,
        )
    else:
        run_pipeline(
            bmrb_id=args.bmrb,
            pdb_id=pdb_id,
            run_dssp=not args.no_dssp,
            save_results=not args.no_save,
        )


def run_batch(bmrb_ids: list, pdb_map: dict = None, run_dssp: bool = True) -> pd.DataFrame:
    """
    Run pipeline on multiple BMRB entries and aggregate statistics.
    pdb_map: {bmrb_id -> pdb_id}
    Also saves per-entry raw merged CSVs to data/batch_cache/ for ML retraining.
    """
    import os
    cache_dir = DATA_DIR / "batch_cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    all_stats = []

    for bmrb_id in bmrb_ids:
        pdb_id = (pdb_map or {}).get(bmrb_id)
        print(f"\n{'#'*60}\n# Processing BMRB {bmrb_id}\n{'#'*60}")
        try:
            res = run_pipeline(bmrb_id, pdb_id=pdb_id, run_dssp=run_dssp, save_results=False)
            if 'stats_df' in res:
                df = res['stats_df'].copy()
                df['bmrb_id'] = bmrb_id
                all_stats.append(df)
            if 'merged_df' in res:
                raw_path = cache_dir / f"bmr{bmrb_id}_raw.csv"
                res['merged_df'].to_csv(raw_path, index=False)
                print(f"  Cached raw shifts → {raw_path}")
        except Exception as e:
            print(f"ERROR on BMRB {bmrb_id}: {e}")
            continue

    if not all_stats:
        return pd.DataFrame()

    combined = pd.concat(all_stats, ignore_index=True)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    out_path = RESULTS_DIR / "combined_stats.csv"
    combined.to_csv(out_path, index=False)
    print(f"\nBatch complete. Combined stats: {len(combined)} rows → {out_path}")
    return combined
