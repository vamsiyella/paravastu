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
    # ── Original confirmed working entries ────────────────────────────────
    (25123, "1UBQ",  "Ubiquitin MPD crystal"),
    (15156, "2LGI",  "GB1 MAS structure"),
    (15283, "2OED",  "GB3 domain"),
    (19025, "1TXQ",  "CAP-Gly 19.9T dataset"),
    (15380, "1PGB",  "GB1 crystal form B1"),
    (17561, "2LBH",  "EETI-II knottin"),
    (11512, "3ONS",  "Ubiquitin alt dataset"),
    (16327, "1FVK",  "DsbA oxidized"),
    (18024, "2K0G",  "CNBD domain"),
    (18397, "1GB1",  "GB1 proton-detected"),
    (19031, "2MPX",  "CAP-Gly+EB1 complex"),
    (25005, "2MPX",  "CAP-Gly on microtubule"),
    # ── Added from scraper v1/v2 ──────────────────────────────────────────
    (5757,  "2AK7",  "Crh-HPr mixed alpha/beta"),
    (17700, "1KEB",  "Thioredoxin"),
    (16964, "2WVN",  "ssNMR helical protein"),
    (50110, "4NUT",  "Snu13p RNA-binding protein"),
    (18808, "2M0G",  "ssNMR mixed alpha/coil"),
    (53330, "9I2I",  "ssNMR mixed strand/coil"),
    (15818, "1HG7",  "Antifreeze protein"),
    (16448, "1PKS",  "BPTI Kunitz domain"),
    (16318, "2JWU",  "Ubiquitin microcrystals"),
    # ── New entries from scraper v3 and literature ────────────────────────
    
    (18493, "2ZUQ",  "ssNMR mixed repeat protein"),
    (18108, "2LME",  "ssNMR helix-rich protein"),
    (25642, "2N3D",  "Bactofilin BacA — pure beta-helix"),
    (6351,  "3ODV",  "ssNMR mixed helix/strand peptide"),
    (16060, "1NOR",  "NOR1 — pure strand/coil"),
    (12019, "1ED7",  "ssNMR pure strand protein"),
    (34178, "6EKA",  "HELLF prion amyloid — new from literature"),
    (25334, "2N7H",  "FimA pilus subunit — Ig-like all-beta"),
    (19747, "2MJZ",  "M13 bacteriophage G8P — helix-rich capsid"),
    (25076, "2MS8",  "MAVS CARD domain — 90% helix"),
    (25788, "2N70",  "Influenza M2 transmembrane helix"),
    (30121, "5KK3",  "Amyloid-beta 42 fibril"),
    (30304, "5W3N",  "FUS LC domain fibril"),
    (18170, "2LNL",  "ssNMR large helix-rich protein"),
    (50411, "7JK8",  "ssNMR recent helix-rich deposit"),
    (30094, "5JZR",  "AP205 phage coat protein"),
]
SKIP_ENTRIES: set = set()


def run_pipeline_batch(
    entries: list = None,
    skip: set = None,
    run_dssp: bool = True,
) -> pd.DataFrame:
    """
    Run the pipeline on all SOLID_STATE_ENTRIES and cache results.
 
    Outputs per-entry files to data/batch_cache/:
        merged_shifts_{bmrb_id}.csv   — shifts + SS labels
        stats_{bmrb_id}.csv           — aggregated statistics
 
    Also writes:
        data/reference_db.csv         — combined stats across all entries
        data/batch_log.csv            — run summary with success/fail per entry
 
    After this completes, run:
        python src/retrainModel.py
    to rebuild the ML models from the full cached dataset.
    """
    
    import time
    if entries is None:
        entries = SOLID_STATE_ENTRIES
    if skip is None:
        skip = SKIP_ENTRIES
 
    batch_cache = DATA_DIR / "batch_cache"
    batch_cache.mkdir(parents=True, exist_ok=True)
 
    all_stats = []
    log_rows = []
 
    print(f"\n{'='*60}")
    print(f"  BATCH RUN — {len(entries)} entries ({len(skip)} skipped)")
    print(f"{'='*60}")
 
    for bmrb_id, pdb_id, label in entries:
        if bmrb_id in skip:
            print(f"\n[SKIP] BMRB {bmrb_id} ({label})")
            continue
 
        print(f"\n{'#'*60}")
        print(f"# BMRB {bmrb_id} / PDB {pdb_id}  — {label}")
        print(f"{'#'*60}")
 
        t0 = time.time()
        try:
            res = run_pipeline(
                bmrb_id=bmrb_id,
                pdb_id=pdb_id,
                run_dssp=run_dssp,
                save_results=False,
            )
 
            if not res or 'merged_df' not in res:
                raise ValueError("Pipeline returned no merged_df")
 
            merged = res['merged_df']
            labeled_frac = (
                merged[merged['ss_class'].isin(['helix','strand','coil'])].shape[0]
                / max(len(merged), 1)
            )
 
            # Cache per-entry files
            merged.to_csv(batch_cache / f"merged_shifts_{bmrb_id}.csv", index=False)
 
            if 'stats_df' in res:
                stats = res['stats_df'].copy()
                stats['bmrb_id'] = bmrb_id
                stats.to_csv(batch_cache / f"stats_{bmrb_id}.csv", index=False)
                all_stats.append(stats)
 
            elapsed = time.time() - t0
            ss_dist = merged['ss_class'].value_counts().to_dict()
            status = "ok" if labeled_frac >= 0.50 else "warn:low_labels"
 
            log_rows.append({
                'bmrb_id': bmrb_id, 'pdb_id': pdb_id, 'label': label,
                'status': status, 'n_shifts': len(merged),
                'labeled_frac': round(labeled_frac, 3),
                'helix': ss_dist.get('helix', 0),
                'strand': ss_dist.get('strand', 0),
                'coil': ss_dist.get('coil', 0),
                'elapsed_s': round(elapsed, 1),
            })
            print(f"  ✓ Done in {elapsed:.0f}s  |  labeled: {labeled_frac:.0%}  |  SS: {ss_dist}")
 
        except Exception as e:
            elapsed = time.time() - t0
            print(f"  ✗ FAILED: {e}")
            log_rows.append({
                'bmrb_id': bmrb_id, 'pdb_id': pdb_id, 'label': label,
                'status': f'error:{type(e).__name__}', 'n_shifts': 0,
                'labeled_frac': 0, 'helix': 0, 'strand': 0, 'coil': 0,
                'elapsed_s': round(elapsed, 1),
            })
 
    # Save batch log
    log_df = pd.DataFrame(log_rows)
    log_df.to_csv(DATA_DIR / "batch_log.csv", index=False)
    print(f"\n[BATCH] Log saved → data/batch_log.csv")
 
    # Save combined reference DB
    if all_stats:
        combined = pd.concat(all_stats, ignore_index=True)
        combined.to_csv(DATA_DIR / "reference_db.csv", index=False)
        n_ok = log_df[log_df['status'] == 'ok'].shape[0]
        n_fail = log_df[log_df['status'].str.startswith('error')].shape[0]
        print(f"[BATCH] Reference DB: {len(combined)} stat rows from {len(all_stats)} entries")
        print(f"[BATCH] {n_ok} succeeded / {n_fail} failed")
        print(f"[BATCH] Saved → data/reference_db.csv")
        print(f"\nNext step: python src/retrainModel.py")
        return combined
 
    print("[BATCH] No stats collected — check errors above.")
    return pd.DataFrame()


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

    



    # Phase  batch args
    batch_group = parser.add_argument_group("Phase 3 — Batch processing")
    batch_group.add_argument(
        "--batch", action="store_true",
        help="Run Phase 3: process curated BMRB entries and build reference DB"
    )
    batch_group.add_argument(
        "--entries", nargs="*", metavar="BMRB:PDB",
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
        custom_entries = None
        if args.entries:
            custom_entries = []
            for pair in args.entries:
                parts = pair.split(":")
                if len(parts) != 2:
                    print(f"Error: --entries format is BMRB_ID:PDB_ID, got '{pair}'")
                    sys.exit(1)
                try:
                    bmrb_int = int(parts[0])
                    pdb_str  = parts[1].upper()
                    match = [(b, p, l) for b, p, l in SOLID_STATE_ENTRIES if b == bmrb_int]
                    label = match[0][2] if match else f"custom {pdb_str}"
                    custom_entries.append((bmrb_int, pdb_str, label))
                except ValueError:
                    print(f"Error: BMRB ID must be an integer, got '{parts[0]}'")
                    sys.exit(1)
        run_pipeline_batch(entries=custom_entries)
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
