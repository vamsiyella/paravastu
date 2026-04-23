"""
pipeline.py — Main orchestrator for the Paravastu NMR structural annotation pipeline.

Usage:
    python pipeline.py                    # run default (BMRB 17561 + PDB 2LBH)
    python pipeline.py --bmrb 17561 --pdb 2LBH
    python pipeline.py --bmrb 17561 --no-dssp   # skip DSSP if not installed
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
from ml_module import run_ml_pipeline


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

    # If no DSSP data but CSV already had real labels, keep them
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
# Main pipeline
# ---------------------------------------------------------------------------

def run_pipeline(
    bmrb_id: int = 17561,
    pdb_id: str = "2LBH",
    run_dssp: bool = True,
    save_results: bool = True,
) -> dict:
    """
    Full pipeline run. Returns dict with all results.
    """
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    results = {}

    # ── Step 1: Load BMRB data ─────────────────────────────────────────────
    print_section(f"STEP 1 — Loading BMRB {bmrb_id}")

    # Try local cache first — validate it's real NMR-STAR, not a stale error page
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
            print("  Or call: run_pipeline_from_csv(csv_path, pdb_id=...)")
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

    # Auto-detect PDB ID if not specified
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

    # ── Step 4: DSSP (optional) ────────────────────────────────────────────
    dssp_map = {}
    seq_mapping = None

    if run_dssp:
        print_section(f"STEP 4 — DSSP from PDB {pdb_id}")

        # Check DSSP available
        mkdssp = find_mkdssp()
        if mkdssp is None:
            print("WARNING: mkdssp not found. Skipping DSSP.")
            print("Fix: conda install -c salilab dssp")
            run_dssp = False

    if run_dssp and pdb_id:
        pdb_path = download_pdb(pdb_id, DATA_DIR)
            
        # Try mkdssp first, fall back to PDB HELIX/SHEET records
        dssp_df = None
        try:
            dssp_df = extract_dssp_full(pdb_path)
            print(f"[DSSP] mkdssp succeeded.")
        except Exception as e:
            print(f"[DSSP] mkdssp failed: {e}")
            print("[DSSP] Falling back to PDB HELIX/SHEET records...")
            try:
                dssp_df = extract_ss_from_pdb_records(pdb_path)
                print("[DSSP] PDB record fallback succeeded.")
                print("NOTE: To get full DSSP data: conda install -c conda-forge libcifpp")
            except Exception as e2:
                print(f"[DSSP] PDB fallback also failed: {e2}")
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
        # Build properly aligned BMRB-indexed SS map using sequence alignment
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

    # ── Step 8: Visualizations (auto-generated every run) ────────────────
    print_section("STEP 8 — Generating plots")
    try:
        results['bmrb_id'] = bmrb_id
        generate_all_plots(results, output_dir=RESULTS_DIR)
    except Exception as e:
        print(f"[VIZ] Plot generation failed: {e}")

    # ── Step 9: ML model training ──────────────────────────────────────────
    if merged_df is not None and not merged_df.empty:
        labeled = merged_df[merged_df['ss_class'].isin(['helix','strand','coil'])]
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
            print(f"  Run batch mode with more BMRB entries to build training data.")

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
# Batch processing (multiple BMRB entries)
# ---------------------------------------------------------------------------

def run_batch(bmrb_ids: list, pdb_map: dict = None, run_dssp: bool = True) -> pd.DataFrame:
    """
    Run pipeline on multiple BMRB entries and aggregate statistics.
    pdb_map: {bmrb_id -> pdb_id} — optional override for PDB lookups.
    Returns combined statistics DataFrame.
    """
    all_stats = []

    for bmrb_id in bmrb_ids:
        pdb_id = (pdb_map or {}).get(bmrb_id)
        print(f"\n{'#'*60}")
        print(f"# Processing BMRB {bmrb_id}")
        print(f"{'#'*60}")
        try:
            res = run_pipeline(bmrb_id, pdb_id=pdb_id, run_dssp=run_dssp, save_results=False)
            if 'stats_df' in res:
                df = res['stats_df'].copy()
                df['bmrb_id'] = bmrb_id
                all_stats.append(df)
        except Exception as e:
            print(f"ERROR on BMRB {bmrb_id}: {e}")
            continue

    if not all_stats:
        return pd.DataFrame()

    combined = pd.concat(all_stats, ignore_index=True)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    combined.to_csv(RESULTS_DIR / "combined_stats.csv", index=False)
    print(f"\nBatch complete. Combined stats: {len(combined)} rows")
    return combined


# ---------------------------------------------------------------------------
# CSV-based entry point (when .str file is unavailable or you already have a CSV)
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

    csv_path: path to CSV with columns: seq_id, residue, atom, shift
              (ss_class column optional — will be added from DSSP if run_dssp=True)

    Use this when:
    - BMRB is unreachable (network issues)
    - You already have your shifts exported to CSV
    - You want to use data from a different source

    Example:
        results = run_pipeline_from_csv(
            csv_path='data/my_shifts.csv',
            pdb_id='2LBH',
            sequence='VLDLDVRT...',
        )
    """
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    results = {'bmrb_id': bmrb_id}

    # ── Load CSV ──────────────────────────────────────────────────────────
    print_section(f"Loading shifts from CSV: {csv_path}")
    shifts_df = pd.read_csv(csv_path)

    # Normalise ss_class codes if present (H→helix, E→strand, C→coil)
    if 'ss_class' in shifts_df.columns:
        coarse = {'H': 'helix', 'G': 'helix', 'I': 'helix',
                  'E': 'strand', 'B': 'strand',
                  'C': 'coil', 'T': 'coil', 'S': 'coil', '-': 'coil',
                  'helix': 'helix', 'strand': 'strand', 'coil': 'coil',
                  'unknown': 'unknown'}
        shifts_df['ss_class'] = shifts_df['ss_class'].map(coarse).fillna('unknown')
    else:
        shifts_df['ss_class'] = 'unknown'

    # Ensure residue one-letter column exists
    if 'residue' not in shifts_df.columns and 'residue_name' in shifts_df.columns:
        THREE_TO_ONE = {
            'ALA':'A','CYS':'C','ASP':'D','GLU':'E','PHE':'F','GLY':'G',
            'HIS':'H','ILE':'I','LYS':'K','LEU':'L','MET':'M','ASN':'N',
            'PRO':'P','GLN':'Q','ARG':'R','SER':'S','THR':'T','VAL':'V',
            'TRP':'W','TYR':'Y',
        }
        shifts_df['residue'] = shifts_df['residue_name'].apply(
            lambda x: THREE_TO_ONE.get(str(x).upper(), x[0] if len(str(x)) == 1 else 'X')
        )

    # ── Detect CSV type ───────────────────────────────────────────────────
    # Two valid CSV formats:
    #   STATS CSV:  columns = residue, atom, ss_class, count, mean, median, std, min, max
    #               (output of the stats step — already aggregated, no seq_id)
    #   RAW CSV:    columns = seq_id, residue, atom, shift, ss_class
    #               (one row per observed shift)
    lbl = bmrb_id or Path(csv_path).stem
    is_stats_csv = 'mean' in shifts_df.columns and 'seq_id' not in shifts_df.columns

    if is_stats_csv:
        print("Detected: STATS CSV (pre-aggregated — mean/std per residue/atom)")
        print(f"  Rows: {len(shifts_df)}, atoms: {sorted(shifts_df['atom'].unique())}")
        print(f"  Residues: {sorted(shifts_df['residue'].unique())}")
        stats_df = shifts_df.copy()
        # Normalise ss_class
        coarse = {'H':'helix','G':'helix','I':'helix','E':'strand','B':'strand',
                  'C':'coil','T':'coil','S':'coil','-':'coil',' ':'coil',
                  'helix':'helix','strand':'strand','coil':'coil','unknown':'unknown'}
        stats_df['ss_class'] = stats_df['ss_class'].map(coarse).fillna('unknown')
        merged_df = None   # no raw rows to merge
        coverage_df = None
    else:
        print("Detected: RAW SHIFTS CSV (one row per observed shift)")
        # Normalise ss_class
        coarse = {'H':'helix','G':'helix','I':'helix','E':'strand','B':'strand',
                  'C':'coil','T':'coil','S':'coil','-':'coil',' ':'coil',
                  'helix':'helix','strand':'strand','coil':'coil','unknown':'unknown'}
        if 'ss_class' in shifts_df.columns:
            shifts_df['ss_class'] = shifts_df['ss_class'].map(coarse).fillna('unknown')
        else:
            shifts_df['ss_class'] = 'unknown'
        print(f"  Rows: {len(shifts_df)}, SS: {shifts_df['ss_class'].unique().tolist()}")
        stats_df = None
        merged_df = None
        coverage_df = None

    results['shifts_df'] = shifts_df
    results['sequence']  = sequence
    results['pdb_id']    = pdb_id

    # ── Coverage (raw CSV only) ────────────────────────────────────────────
    if not is_stats_csv:
        print_section("Coverage analysis")
        coverage_df = compute_coverage(shifts_df)
        print(f"Residues with data: {len(coverage_df)}, avg atoms/res: {coverage_df['n_atoms'].mean():.2f}")
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

            # Try mkdssp first, fall back to PDB HELIX/SHEET records
            dssp_df = None
            try:
                dssp_df = extract_dssp_full(pdb_path_dl)
                print(f"[DSSP] mkdssp succeeded.")
            except Exception as e:
                print(f"[DSSP] mkdssp failed: {e}")
                print("[DSSP] Falling back to PDB HELIX/SHEET records...")
                try:
                    dssp_df = extract_ss_from_pdb_records(pdb_path_dl)
                    print("[DSSP] PDB record fallback succeeded.")
                    print("NOTE: PDB fallback gives helix/strand/coil only (no phi/psi/accessibility).")
                    print("      For full DSSP data, fix mkdssp with:")
                    print("        conda install -c conda-forge libcifpp")
                    print("      then re-run the pipeline.")
                except Exception as e2:
                    print(f"[DSSP] PDB fallback also failed: {e2}")

            if dssp_df is not None and not dssp_df.empty:
                dssp_map    = dict(zip(dssp_df['residue_number'], dssp_df['ss_class']))
                segments    = build_ss_segments(dssp_map)
                results['dssp_df']  = dssp_df
                results['segments'] = segments
                print(f"SS labels: {len(dssp_map)} residues")
                print(dssp_df['ss_class'].value_counts().to_string())
                for seg in segments:
                    print(f"  {seg['ss_class']:7s} {seg['start']:3d}-{seg['end']:3d} ({seg['length']} aa)")
                if sequence:
                    pdb_seq     = get_pdb_sequence_from_file(pdb_path_dl)
                    seq_mapping = align_and_map(sequence, pdb_seq)
                if save_results:
                    dssp_df.to_csv(RESULTS_DIR / f"dssp_{pdb_id}.csv", index=False)

    # ── Merge (raw CSV) / annotate stats (stats CSV) ──────────────────────
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
        # Stats CSV: if DSSP ran, we can annotate which residues belong to which SS
        # and tell the user to re-run from BMRB for full integration.
        # For now use the stats as-is for prediction.
        print_section("Using pre-computed statistics")
        stats_df = add_random_coil_deviation(stats_df)
        results['stats_df'] = stats_df
        if save_results:
            stats_df.to_csv(RESULTS_DIR / f"stats_{lbl}.csv", index=False)
        print(f"{len(stats_df)} stat rows (all ss_class=unknown)")
        if dssp_map:
            print("NOTE: DSSP ran successfully. To get SS-aware stats, run:")
            print(f"  python src/pipeline.py --bmrb {bmrb_id or 17561} --pdb {pdb_id}")
            print("  This will fetch raw shifts from BMRB and merge with DSSP labels.")

    # ── Prediction ────────────────────────────────────────────────────────
    if sequence and stats_df is not None and not stats_df.empty:
        print_section("Shift prediction")
        predictor = ShiftPredictor(stats_df)
        # Build a properly aligned BMRB-indexed SS map
        # (dssp_map keys are DSSP sequential numbers, not BMRB seq_ids)
        if dssp_map and results.get('dssp_df') is not None and seq_mapping:
            ss_map_for_pred = build_prediction_ss_map(
                results['dssp_df'], seq_mapping, len(sequence)
            )
        elif dssp_map:
            ss_map_for_pred = {i+1: dssp_map.get(i+1, 'coil') for i in range(len(sequence))}
        else:
            ss_map_for_pred = {i+1: 'unknown' for i in range(len(sequence))}
        preds = predictor.predict(sequence, ss_map_for_pred, atoms=['CA', 'N'])
        preds_dev = preds.rename(columns={'predicted_shift': 'shift'})
        preds_dev = add_random_coil_deviation(preds_dev)
        preds_dev.rename(columns={'shift': 'predicted_shift'}, inplace=True)
        results['predictions'] = preds_dev
        if save_results:
            preds_dev.to_csv(RESULTS_DIR / f"predictions_{lbl}.csv", index=False)
        print(f"Predicted shifts for {len(preds_dev)} (residue, atom) pairs")
        print(preds_dev.head(10).to_string())

    # ── Visualizations ────────────────────────────────────────────────────
    print_section("Generating plots")
    try:
        results['bmrb_id'] = lbl
        generate_all_plots(results, output_dir=RESULTS_DIR)
    except Exception as e:
        print(f"[VIZ] Plot generation failed: {e}")

    # ── ML (raw CSV only — needs per-shift rows) ───────────────────────────
    if not is_stats_csv and merged_df is not None:
        labeled = merged_df[merged_df['ss_class'].isin(['helix','strand','coil'])]
        if len(labeled) >= 15 and labeled['ss_class'].nunique() >= 2:
            print_section("Training ML models")
            try:
                ml_results = run_ml_pipeline(merged_df, results_dir=RESULTS_DIR)
                results['ml'] = ml_results
            except Exception as e:
                print(f"[ML] Training failed: {e}")

    # ── Summary ───────────────────────────────────────────────────────────
    print_section("PIPELINE COMPLETE")
    print(f"  Source CSV:     {csv_path}")
    print(f"  CSV type:       {'stats (pre-aggregated)' if is_stats_csv else 'raw shifts'}")
    print(f"  PDB:            {pdb_id or 'none'}")
    print(f"  DSSP residues:  {len(dssp_map)}")
    print(f"  Stat rows:      {len(stats_df) if stats_df is not None else 0}")
    print(f"  Plots:    {len(list(RESULTS_DIR.glob('*.png')))} PNG files")
    print(f"  Results → {RESULTS_DIR}")

    return results


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Paravastu NMR Structural Annotation Pipeline"
    )
    parser.add_argument("--bmrb",     type=int,  default=17561, help="BMRB entry ID")
    parser.add_argument("--pdb",      type=str,  default="2LBH", help="PDB ID (or 'auto')")
    parser.add_argument("--csv",      type=str,  default=None,  help="Path to shifts CSV (bypasses BMRB fetch)")
    parser.add_argument("--seq",      type=str,  default=None,  help="Protein sequence (one-letter, used with --csv)")
    parser.add_argument("--no-dssp",  action="store_true",      help="Skip DSSP extraction")
    parser.add_argument("--no-save",  action="store_true",      help="Don't save result files")
    args = parser.parse_args()

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
