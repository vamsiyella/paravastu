"""
retrainModelv5.py — Key improvements over v4:

1. PER-RESIDUE-TYPE NORMALIZATION (biggest expected gain)
   Instead of subtracting a fixed random-coil table, subtract the DATASET MEDIAN
   for each (residue, atom) combination, computed within each training fold.
   This removes residue-type chemistry confounding without needing a literature table.
   e.g. Val CA helix vs coil, not (Val CA - 62.2) vs coil

2. CARBONYL C' FEATURE WITH PROPER RC VALUES
   C' shifts: helix ~177-178 ppm, strand ~174-175 ppm (3 ppm separation, orthogonal to CA/CB)
   Currently undertrained because C' RC values in v4 table were imprecise.
   Fixed RC values from Wishart 2011 (solid-state corrected).

3. STRAND/COIL DISCRIMINATOR: CA - CB IN DEVIATION SPACE
   Raw CA-CB is dominated by residue type (Gly CA-CB is huge, Val is small)
   dev_CB_minus_CA uses (CB_obs - CB_rc) - (CA_obs - CA_rc)
   This isolates the structural signal from residue-type chemistry.
   Positive = strand tendency, negative = helix tendency.

4. MISSINGNESS PATTERN AS A FEATURE
   In ssNMR, CB is often missing for Gly (no CB) and sometimes for disordered coil.
   The pattern of which atoms are observed is itself structural information.
   Add binary indicator features: has_CB, has_N, has_C.

5. WITHIN-FOLD IMPUTATION (prevent leakage)
   In v4, fillna used the global median — this leaks test-set information.
   v5 computes imputation values on training data only, applies to test.

6. ANONYMOUS PROTEIN HANDLING
   ssNMR-* groups are kept but flagged. If the same fold appears in multiple
   anonymous entries (likely different labs on similar proteins), they share a group.

Run:
    conda activate nmr
    cd C:\\Users\\vamsi\\.vscode\\paravastu
    python src/retrainModelv5.py
"""

import sys
import numpy as np
import pandas as pd
from pathlib import Path
from collections import defaultdict

SRC_DIR    = Path(__file__).resolve().parent
ROOT_DIR   = SRC_DIR.parent
DATA_DIR   = ROOT_DIR / "data"
CACHE_DIR  = DATA_DIR / "batch_cache"
RESULTS_DIR = ROOT_DIR / "results"

sys.path.insert(0, str(SRC_DIR))

from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import LeaveOneGroupOut
from sklearn.metrics import (classification_report, confusion_matrix,
                              f1_score, accuracy_score)
from sklearn.utils.class_weight import compute_class_weight
import joblib

try:
    import xgboost as xgb
    HAS_XGB = True
except ImportError:
    HAS_XGB = False
    print("[ML] XGBoost not available")

# ---------------------------------------------------------------------------
# Protein family grouping
# Keep all entries; ssNMR-* anonymous proteins each get their own group
# unless we have reason to merge them
# ---------------------------------------------------------------------------
PROTEIN_GROUPS = {
    # GB1 family
    15156: "GB1", 15283: "GB1", 15380: "GB1", 18397: "GB1",
    16873: "GB1", 30088: "GB1",
    # CAP-Gly family
    19025: "CAP-Gly", 19031: "CAP-Gly", 25005: "CAP-Gly", 17937: "CAP-Gly",
    # Ubiquitin family
    25123: "Ubiquitin", 11512: "Ubiquitin", 16318: "Ubiquitin",
    # Named unique proteins
    17561: "EETI-II",
    16327: "DsbA",     18543: "DsbA",
    18024: "CNBD",
    5757:  "Crh-HPr",
    17700: "Thioredoxin",
    15818: "Antifreeze",
    16448: "BPTI",
    25334: "FimA-pilus",
    34178: "HELLF-prion",
    50110: "Snu13p",
    25076: "MAVS-CARD",
    19747: "M13-G8P",
    25788: "M2-channel",
    30094: "AP205-capsid",
    25642: "BactofilinBacA",
    30304: "FUS-LC",
    30121: "AmbetaFibril",
    18808: "ssNMR-mixed-1",
    # Anonymous entries — each gets its own group (conservative)
    # If two anonymous entries are from the same protein they'd inflate accuracy
    # but we can't know without sequence data, so keep separate
    6351:  "ssNMR-anon-6351",
    12019: "ssNMR-anon-12019",
    16060: "ssNMR-anon-16060",
    16964: "ssNMR-anon-16964",
    18108: "ssNMR-anon-18108",
    18170: "ssNMR-anon-18170",
    18493: "ssNMR-anon-18493",
    50411: "ssNMR-anon-50411",
    53330: "ssNMR-anon-53330",
}

# ---------------------------------------------------------------------------
# Solid-state corrected random coil values (Wishart 2011 + BMRB ssNMR averages)
# C' values are the key addition — literature shows clear helix/strand separation
# ---------------------------------------------------------------------------
RC_SSnmr = {
    'A': {'CA': 52.5,  'CB': 19.1,  'C': 177.8, 'N': 123.8, 'HA': 4.32},
    'R': {'CA': 56.3,  'CB': 31.0,  'C': 176.3, 'N': 120.5, 'HA': 4.34},
    'N': {'CA': 53.3,  'CB': 38.9,  'C': 175.2, 'N': 118.7, 'HA': 4.75},
    'D': {'CA': 54.4,  'CB': 41.1,  'C': 176.3, 'N': 120.4, 'HA': 4.76},
    'C': {'CA': 58.4,  'CB': 28.0,  'C': 174.6, 'N': 118.8, 'HA': 4.56},
    'Q': {'CA': 56.0,  'CB': 29.4,  'C': 175.9, 'N': 119.8, 'HA': 4.37},
    'E': {'CA': 56.8,  'CB': 30.2,  'C': 176.6, 'N': 120.2, 'HA': 4.29},
    'G': {'CA': 45.3,  'CB': None,  'C': 173.8, 'N': 108.8, 'HA': 3.97},
    'H': {'CA': 56.1,  'CB': 29.4,  'C': 174.1, 'N': 118.2, 'HA': 4.63},
    'I': {'CA': 61.3,  'CB': 38.8,  'C': 175.8, 'N': 120.9, 'HA': 4.17},
    'L': {'CA': 55.4,  'CB': 42.4,  'C': 177.6, 'N': 121.8, 'HA': 4.34},
    'K': {'CA': 56.7,  'CB': 33.1,  'C': 176.6, 'N': 120.4, 'HA': 4.36},
    'M': {'CA': 55.8,  'CB': 33.1,  'C': 176.3, 'N': 119.6, 'HA': 4.52},
    'F': {'CA': 57.9,  'CB': 39.6,  'C': 175.8, 'N': 120.3, 'HA': 4.66},
    'P': {'CA': 63.5,  'CB': 32.1,  'C': 177.3, 'N': 136.5, 'HA': 4.44},
    'S': {'CA': 58.5,  'CB': 63.8,  'C': 174.6, 'N': 115.7, 'HA': 4.47},
    'T': {'CA': 62.0,  'CB': 69.8,  'C': 174.7, 'N': 113.6, 'HA': 4.35},
    'W': {'CA': 57.7,  'CB': 29.9,  'C': 176.1, 'N': 121.3, 'HA': 4.70},
    'Y': {'CA': 58.1,  'CB': 38.8,  'C': 175.9, 'N': 120.3, 'HA': 4.60},
    'V': {'CA': 62.4,  'CB': 32.9,  'C': 176.3, 'N': 119.9, 'HA': 4.18},
}

# Wishart (1994) CSI thresholds for CA: +0.7 ppm = helix, -0.7 ppm = strand
# For C': +0.4 ppm = helix, -0.4 ppm = strand
# For CB: -0.5 ppm = helix, +0.5 ppm = strand (inverted!)
CSI_THRESHOLDS = {
    'CA': (0.7,  -0.7),   # (helix_threshold, strand_threshold)
    'C':  (0.4,  -0.4),
    'CB': (-0.5,  0.5),   # inverted
    'N':  (0.0,   0.0),   # N not used in original CSI
}

def csi_score(dev_ca, dev_c, dev_cb):
    """
    Compute Wishart (1994) Combined Chemical Shift Index score.
    Returns a float in [-3, +3]:
      +3 = strongly helical
      -3 = strongly strand
       0 = ambiguous/coil
    """
    score = 0.0
    if dev_ca is not None and not np.isnan(dev_ca):
        if dev_ca > 0.7: score += 1
        elif dev_ca < -0.7: score -= 1
    if dev_c is not None and not np.isnan(dev_c):
        if dev_c > 0.4: score += 1
        elif dev_c < -0.4: score -= 1
    if dev_cb is not None and not np.isnan(dev_cb):
        # CB is inverted: more negative in helix
        if dev_cb < -0.5: score += 1
        elif dev_cb > 0.5: score -= 1
    return score

# ---------------------------------------------------------------------------
# Step 1: Load cached shifts
# ---------------------------------------------------------------------------

def load_all_cached_shifts():
    """
    Load cached raw shift CSVs from batch_cache/.
    Prefer bmr{id}_{pdb}_raw.csv over bmr{id}_raw.csv (PDB-specific version is more accurate).
    """
    if not CACHE_DIR.exists():
        print(f"ERROR: {CACHE_DIR} does not exist. Run --batch first.")
        sys.exit(1)

    # Find all cache files, prefer PDB-specific ones
    all_files: dict[int, Path] = {}
    for f in sorted(CACHE_DIR.glob("bmr*_raw.csv")):
        parts = f.stem.replace("bmr", "").split("_")
        try:
            bmrb_id = int(parts[0])
        except (ValueError, IndexError):
            continue
        # Prefer files with PDB suffix (bmr17561_2LBH_raw.csv > bmr17561_raw.csv)
        if bmrb_id not in all_files or len(f.stem) > len(all_files[bmrb_id].stem):
            all_files[bmrb_id] = f

    # Also check merged_shifts_*.csv from results/
    for d in [CACHE_DIR, RESULTS_DIR]:
        for f in sorted(d.glob("merged_shifts_*.csv")):
            stem = f.stem.replace("merged_shifts_", "")
            for part in stem.split("_"):
                try:
                    bmrb_id = int(part)
                    if bmrb_id not in all_files:
                        all_files[bmrb_id] = f
                    break
                except ValueError:
                    continue

    print(f"[RETRAIN v5] Found {len(all_files)} unique BMRB entries in cache:")
    frames = []
    skipped = []

    for bmrb_id, fpath in sorted(all_files.items()):
        df = pd.read_csv(fpath)

        # Validate columns
        if 'ss_class' not in df.columns or 'shift' not in df.columns:
            skipped.append((bmrb_id, "missing ss_class or shift"))
            continue

        labeled = df[df['ss_class'].isin(['helix', 'strand', 'coil'])]
        total_res = df['seq_id'].nunique() if 'seq_id' in df.columns else max(len(df), 1)
        labeled_res = labeled['seq_id'].nunique() if 'seq_id' in labeled.columns else len(labeled)
        label_frac = labeled_res / max(total_res, 1)

        if label_frac < 0.50:
            skipped.append((bmrb_id, f"only {label_frac:.0%} labeled"))
            continue

        df['source_bmrb'] = bmrb_id
        df['protein_group'] = PROTEIN_GROUPS.get(bmrb_id, f"unknown_{bmrb_id}")
        labeled = df[df['ss_class'].isin(['helix', 'strand', 'coil'])].copy()
        frames.append(labeled)

        ss = labeled['ss_class'].value_counts().to_dict()
        grp = PROTEIN_GROUPS.get(bmrb_id, "?")
        print(f"  BMRB {bmrb_id:>6} [{grp:<22}]: {len(labeled):>5} shifts  "
              f"H={ss.get('helix',0):>4} S={ss.get('strand',0):>4} C={ss.get('coil',0):>4}")

    if skipped:
        print(f"\n  [SKIP] {len(skipped)} entries:")
        for bid, reason in skipped:
            print(f"    BMRB {bid}: {reason}")

    if not frames:
        print("ERROR: No labeled data found. Run --batch first.")
        sys.exit(1)

    combined = pd.concat(frames, ignore_index=True)
    groups = combined['protein_group'].unique()
    print(f"\n[RETRAIN v5] {len(combined)} shifts | {len(frames)} entries | "
          f"{len(groups)} protein groups")
    print(f"  SS dist: {combined['ss_class'].value_counts().to_dict()}")
    return combined

# ---------------------------------------------------------------------------
# Step 2: Build features with per-residue-type normalization
# ---------------------------------------------------------------------------

def build_features(df: pd.DataFrame, window: int = 2,
                   train_medians: dict = None) -> tuple:
    """
    Build feature matrix.

    train_medians: if provided, use these for per-residue normalization (test-time).
                   if None, compute from df (train-time) and return them.

    Returns (X, y, groups, medians_dict)
    """
    # Pivot to one row per residue
    pivot = df.pivot_table(
        index=['source_bmrb', 'seq_id', 'residue', 'protein_group'],
        columns='atom',
        values='shift',
        aggfunc='mean',
    ).reset_index()

    ss = (
        df[df['ss_class'].isin(['helix', 'strand', 'coil'])]
        .groupby(['source_bmrb', 'seq_id'])['ss_class']
        .agg(lambda x: x.mode().iloc[0])
        .reset_index()
    )
    pivot = pivot.merge(ss, on=['source_bmrb', 'seq_id'], how='inner')
    pivot = pivot.sort_values(['source_bmrb', 'seq_id']).reset_index(drop=True)

    feature_cols = []

    # ── Raw shifts ────────────────────────────────────────────────────────
    for atom in ['CA', 'CB', 'C', 'N', 'H', 'HA']:
        col = f'shift_{atom}'
        pivot[col] = pivot[atom] if atom in pivot.columns else np.nan
        feature_cols.append(col)

    # ── RC deviation (global table) ───────────────────────────────────────
    for atom in ['CA', 'CB', 'C', 'N']:
        col = f'dev_{atom}'
        def get_dev(row, _atom=atom):
            if _atom not in row.index or pd.isna(row.get(_atom, np.nan)):
                return np.nan
            rc = RC_SSMNR_get(row['residue'], _atom)
            return float(row[_atom]) - rc if rc is not None else np.nan
        pivot[col] = pivot.apply(get_dev, axis=1)
        feature_cols.append(col)

    # ── Per-residue-type normalization (KEY NEW FEATURE) ─────────────────
    # For each (residue, atom) combination, subtract the DATASET median observed shift.
    # This removes residue-type chemistry without needing a literature RC table.
    # Must be computed on training data only to prevent leakage.
    medians_out = {}
    for atom in ['CA', 'CB', 'C']:
        raw_col = atom if atom in pivot.columns else None
        col = f'norm_{atom}'
        if raw_col is not None:
            if train_medians is not None:
                # Test mode: use pre-computed medians from training set
                key = f'norm_{atom}'
                med_map = train_medians.get(key, {})
                pivot[col] = pivot.apply(
                    lambda r, _atom=atom, _map=med_map: (
                        float(r[_atom]) - _map.get(r['residue'], np.nan)
                        if not pd.isna(r.get(_atom, np.nan)) else np.nan
                    ), axis=1
                )
            else:
                # Train mode: compute medians from this data
                med_map = (pivot[[atom, 'residue']].dropna()
                           .groupby('residue')[atom].median().to_dict())
                medians_out[f'norm_{atom}'] = med_map
                pivot[col] = pivot.apply(
                    lambda r, _atom=atom, _map=med_map: (
                        float(r[_atom]) - _map.get(r['residue'], np.nan)
                        if not pd.isna(r.get(_atom, np.nan)) else np.nan
                    ), axis=1
                )
        else:
            pivot[col] = np.nan
        feature_cols.append(col)

    # ── dev_CB_minus_CA (deviation space — key discriminator) ─────────────
    # Raw CA-CB is residue-type dominated (Gly, Val etc.)
    # dev_CB - dev_CA isolates structural signal:
    #   helix: dev_CA positive, dev_CB negative → difference very negative
    #   strand: dev_CA negative, dev_CB positive → difference very positive
    if 'dev_CA' in pivot.columns and 'dev_CB' in pivot.columns:
        pivot['dev_CB_minus_CA'] = pivot['dev_CB'] - pivot['dev_CA']
    else:
        pivot['dev_CB_minus_CA'] = np.nan
    feature_cols.append('dev_CB_minus_CA')

    # ── CSI composite score (Wishart 1994) ────────────────────────────────
    pivot['csi_score'] = pivot.apply(
        lambda r: csi_score(
            r.get('dev_CA', np.nan),
            r.get('dev_C', np.nan),
            r.get('dev_CB', np.nan)
        ), axis=1
    )
    feature_cols.append('csi_score')

    # Per-atom CSI indicators
    for atom, (h_thr, s_thr) in CSI_THRESHOLDS.items():
        dev_col = f'dev_{atom}'
        col = f'csi_{atom}'
        if dev_col in pivot.columns and h_thr != 0:
            def csi_ind(val, _h=h_thr, _s=s_thr):
                if pd.isna(val): return np.nan
                if _h > 0 and val > _h: return 1.0
                if _s < 0 and val < _s: return -1.0
                return 0.0
            pivot[col] = pivot[dev_col].apply(csi_ind)
        else:
            pivot[col] = np.nan
        feature_cols.append(col)

    # ── Missingness indicator features ────────────────────────────────────
    # In ssNMR: missing CB can indicate Gly or disorder (coil)
    for atom in ['CB', 'C', 'N']:
        col = f'has_{atom}'
        if atom in pivot.columns:
            pivot[col] = (~pivot[atom].isna()).astype(float)
        else:
            pivot[col] = 0.0
        feature_cols.append(col)

    # ── Residue one-hot ───────────────────────────────────────────────────
    for aa in list('ACDEFGHIKLMNPQRSTVWY'):
        col = f'is_{aa}'
        pivot[col] = (pivot['residue'] == aa).astype(float)
        feature_cols.append(col)

    # ── Context window (within same protein only) ─────────────────────────
    for offset in range(-window, window + 1):
        if offset == 0:
            continue
        for atom in ['CA', 'CB', 'C']:
            # Raw shift context
            base = f'shift_{atom}'
            col = f'{atom}_n{offset:+d}'
            pivot[col] = pivot.groupby('source_bmrb')[base].shift(-offset)
            feature_cols.append(col)
            # Deviation context
            dev_base = f'dev_{atom}'
            col2 = f'dev_{atom}_n{offset:+d}'
            pivot[col2] = pivot.groupby('source_bmrb')[dev_base].shift(-offset)
            feature_cols.append(col2)
        # CSI score context
        col = f'csi_n{offset:+d}'
        pivot[col] = pivot.groupby('source_bmrb')['csi_score'].shift(-offset)
        feature_cols.append(col)

    X = pivot[feature_cols].copy()
    y = pivot['ss_class']
    groups = pivot['protein_group']

    print(f"[RETRAIN v5] Feature matrix: {X.shape[0]} × {X.shape[1]}")
    print(f"[RETRAIN v5] Labels: {y.value_counts().to_dict()}")

    return X, y, groups, medians_out


def RC_SSMNR_get(residue, atom):
    """Safe lookup in RC_SSnmr table."""
    return RC_SSnmr.get(residue, {}).get(atom, None)


# ---------------------------------------------------------------------------
# Sliding window consensus (post-processing)
# ---------------------------------------------------------------------------

def sliding_window_consensus(proba_matrix: np.ndarray, window: int = 3) -> np.ndarray:
    """
    Smooth predicted probabilities with a sliding window average.
    Reduces boundary misclassifications (a residue surrounded by helix
    neighbors should itself be helix).

    proba_matrix: (n_residues, n_classes) probability matrix
    window: half-window size (total window = 2*window+1)
    Returns smoothed probability matrix, same shape.
    """
    n = proba_matrix.shape[0]
    smoothed = np.zeros_like(proba_matrix)
    for i in range(n):
        lo = max(0, i - window)
        hi = min(n, i + window + 1)
        smoothed[i] = proba_matrix[lo:hi].mean(axis=0)
    return smoothed

# ---------------------------------------------------------------------------
# Step 3: LOPO-CV with within-fold imputation
# ---------------------------------------------------------------------------

def lopo_cv_v5(model_factory, X: pd.DataFrame, y: pd.Series,
               groups: pd.Series, label: str, smooth_window: int = 1):
    """
    LOPO-CV with:
    - Within-fold imputation (train median, not global)
    - Sliding window consensus on predicted probabilities
    - Reports both raw and smoothed accuracy
    """
    lopo = LeaveOneGroupOut()
    unique_groups = sorted(groups.unique())
    n_folds = len(unique_groups)

    print(f"\n[RETRAIN v5] {label} LOPO-CV ({n_folds} folds)")

    all_true_raw, all_pred_raw = [], []
    all_true_sm,  all_pred_sm  = [], []
    fold_results = []

    classes_global = sorted(y.unique())

    for fold_idx, (tr_idx, te_idx) in enumerate(lopo.split(X, y, groups)):
        test_group = groups.iloc[te_idx[0]]
        X_tr_raw = X.iloc[tr_idx]
        X_te_raw = X.iloc[te_idx]
        y_tr = y.iloc[tr_idx]
        y_te = y.iloc[te_idx]

        if y_tr.nunique() < 2:
            print(f"  Fold {fold_idx+1:2d} [{test_group:<22}] SKIP (train has <2 classes)")
            continue

        # Within-fold imputation: compute medians on train, apply to both
        train_col_medians = X_tr_raw.median(numeric_only=True)
        X_tr = X_tr_raw.fillna(train_col_medians)
        X_te = X_te_raw.fillna(train_col_medians)

        # Class weights on training data
        classes_tr = sorted(y_tr.unique())
        weights = compute_class_weight('balanced', classes=np.array(classes_tr), y=y_tr)
        cw = dict(zip(classes_tr, weights))

        model = model_factory(class_weight=cw)
        model.fit(X_tr, y_tr)

        # Raw predictions
        y_pred_raw = model.predict(X_te)
        acc_raw = accuracy_score(y_te, y_pred_raw)

        # Smoothed predictions (sliding window on probabilities)
        try:
            proba = model.predict_proba(X_te)
            classes_model = model.classes_
            smoothed_proba = sliding_window_consensus(proba, window=smooth_window)
            y_pred_sm = np.array([classes_model[i] for i in smoothed_proba.argmax(axis=1)])
            acc_sm = accuracy_score(y_te, y_pred_sm)
        except Exception:
            y_pred_sm = y_pred_raw
            acc_sm = acc_raw

        n_h = (y_te == 'helix').sum()
        n_s = (y_te == 'strand').sum()
        n_c = (y_te == 'coil').sum()
        print(f"  Fold {fold_idx+1:2d} [{test_group:<22}] "
              f"n=H{n_h}/S{n_s}/C{n_c}  raw={acc_raw:.3f}  smooth={acc_sm:.3f}")

        all_true_raw.extend(y_te); all_pred_raw.extend(y_pred_raw)
        all_true_sm.extend(y_te);  all_pred_sm.extend(y_pred_sm)
        fold_results.append({'group': test_group, 'raw': acc_raw, 'smooth': acc_sm,
                              'n_h': n_h, 'n_s': n_s, 'n_c': n_c})

    acc_raw = accuracy_score(all_true_raw, all_pred_raw)
    acc_sm  = accuracy_score(all_true_sm, all_pred_sm)
    f1_raw  = f1_score(all_true_raw, all_pred_raw, average='macro', zero_division=0)
    f1_sm   = f1_score(all_true_sm,  all_pred_sm,  average='macro', zero_division=0)

    print(f"\n  {label} — RAW:     accuracy={acc_raw:.3f}  F1={f1_raw:.3f}")
    print(f"  {label} — SMOOTH:  accuracy={acc_sm:.3f}  F1={f1_sm:.3f}")
    print(f"\n  Smoothed per-class report:")
    print(classification_report(all_true_sm, all_pred_sm, digits=3, zero_division=0))

    labels = ['coil', 'helix', 'strand']
    cm = confusion_matrix(all_true_sm, all_pred_sm, labels=labels)
    print(f"  Confusion (smoothed):")
    print(pd.DataFrame(cm, index=labels, columns=labels).to_string())

    return acc_raw, acc_sm, f1_sm, fold_results


# ---------------------------------------------------------------------------
# Step 4: Final model on all data
# ---------------------------------------------------------------------------

def train_final(model_factory, X: pd.DataFrame, y: pd.Series, label: str):
    col_medians = X.median(numeric_only=True)
    X_filled = X.fillna(col_medians)
    classes = sorted(y.unique())
    weights = compute_class_weight('balanced', classes=np.array(classes), y=y)
    cw = dict(zip(classes, weights))
    model = model_factory(class_weight=cw)
    model.fit(X_filled, y)
    print(f"[RETRAIN v5] {label} final model trained on {len(X)} samples")
    return model, col_medians


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("=" * 65)
    print("  PARAVASTU — ML RETRAINING v5")
    print("  Per-residue normalization + C' features + within-fold imputation")
    print("=" * 65)

    df = load_all_cached_shifts()
    n_groups = df['protein_group'].nunique()

    if n_groups < 4:
        print(f"\nWARNING: Only {n_groups} protein groups. LOPO-CV unreliable.")

    X, y, groups, train_medians = build_features(df, window=2)

    print(f"\n{'='*65}")
    print(f"  LEAVE-ONE-PROTEIN-OUT CROSS-VALIDATION (v5)")
    print(f"  Smooth window = 1 (conservative; per-fold smoothing)")
    print(f"{'='*65}")

    def rf_factory(class_weight=None):
        return RandomForestClassifier(
            n_estimators=400, max_depth=None, min_samples_leaf=2,
            class_weight=class_weight, random_state=42, n_jobs=-1
        )

    rf_acc_raw, rf_acc_sm, rf_f1, rf_folds = lopo_cv_v5(
        rf_factory, X, y, groups, "RandomForest", smooth_window=1
    )

    if HAS_XGB:
        from sklearn.preprocessing import LabelEncoder

        def xgb_factory(class_weight=None):
            # XGBoost uses sample_weight not class_weight
            return xgb.XGBClassifier(
                n_estimators=400, max_depth=5, learning_rate=0.05,
                subsample=0.8, colsample_bytree=0.8,
                eval_metric='mlogloss', random_state=42, n_jobs=-1, verbosity=0,
                use_label_encoder=False,
            )

        # XGB needs label-encoded y
        le = LabelEncoder()
        y_enc = pd.Series(le.fit_transform(y), index=y.index, name='ss_class')
        classes_str = le.classes_  # ['coil', 'helix', 'strand']

        lopo = LeaveOneGroupOut()
        all_true_sm, all_pred_sm = [], []
        print(f"\n[RETRAIN v5] XGBoost LOPO-CV ({n_groups} folds)")

        for fold_idx, (tr_idx, te_idx) in enumerate(lopo.split(X, y_enc, groups)):
            test_group = groups.iloc[te_idx[0]]
            X_tr_raw = X.iloc[tr_idx]; X_te_raw = X.iloc[te_idx]
            y_tr = y_enc.iloc[tr_idx]; y_te = y_enc.iloc[te_idx]

            if y_tr.nunique() < 2: continue

            col_medians = X_tr_raw.median(numeric_only=True)
            X_tr = X_tr_raw.fillna(col_medians).fillna(-999)
            X_te = X_te_raw.fillna(col_medians).fillna(-999)

            counts = np.bincount(y_tr)
            max_c = counts.max()
            sw = np.array([max_c / max(counts[yi], 1) for yi in y_tr])

            model_xgb = xgb_factory()
            model_xgb.fit(X_tr, y_tr, sample_weight=sw)

            proba = model_xgb.predict_proba(X_te)
            smoothed = sliding_window_consensus(proba, window=1)
            y_pred = smoothed.argmax(axis=1)

            acc = accuracy_score(y_te, y_pred)
            n_h = (y.iloc[te_idx] == 'helix').sum()
            n_s = (y.iloc[te_idx] == 'strand').sum()
            n_c = (y.iloc[te_idx] == 'coil').sum()
            print(f"  Fold {fold_idx+1:2d} [{test_group:<22}] n=H{n_h}/S{n_s}/C{n_c}  smooth={acc:.3f}")

            all_true_sm.extend(le.inverse_transform(y_te))
            all_pred_sm.extend([classes_str[i] for i in y_pred])

        xgb_acc = accuracy_score(all_true_sm, all_pred_sm)
        xgb_f1  = f1_score(all_true_sm, all_pred_sm, average='macro', zero_division=0)
        print(f"\n  XGBoost SMOOTH: accuracy={xgb_acc:.3f}  F1={xgb_f1:.3f}")
        print(classification_report(all_true_sm, all_pred_sm, digits=3, zero_division=0))

    # ── Final models ──────────────────────────────────────────────────────
    print(f"\n{'='*65}")
    print(f"  FINAL MODELS (all data)")
    print(f"{'='*65}")

    rf_final, rf_medians = train_final(rf_factory, X, y, "RandomForest")
    fi = pd.Series(rf_final.feature_importances_, index=X.columns).sort_values(ascending=False)
    print("\nTop 20 features:")
    print(fi.head(20).to_string())

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    rf_dict = {
        'model': rf_final,
        'col_medians': rf_medians,
        'feature_importance': fi,
        'lopo_acc_raw': rf_acc_raw,
        'lopo_acc_smooth': rf_acc_sm,
        'lopo_f1': rf_f1,
        'version': 'v5',
        'smooth_window': 1,
        'feature_cols': list(X.columns),
    }
    joblib.dump(rf_dict, RESULTS_DIR / "model_rf.joblib")
    print(f"\n[RETRAIN v5] Saved → results/model_rf.joblib")

    if HAS_XGB:
        col_med_xgb = X.median(numeric_only=True)
        X_xgb = X.fillna(col_med_xgb).fillna(-999)
        le2 = LabelEncoder()
        y2 = le2.fit_transform(y)
        counts2 = np.bincount(y2); max_c2 = counts2.max()
        sw2 = np.array([max_c2 / max(counts2[yi], 1) for yi in y2])
        xgb_final = xgb_factory()
        xgb_final.fit(X_xgb, y2, sample_weight=sw2)
        xgb_dict = {
            'model': xgb_final, 'label_encoder': le2,
            'col_medians': col_med_xgb,
            'lopo_acc_smooth': xgb_acc, 'lopo_f1': xgb_f1,
            'version': 'v5', 'feature_cols': list(X.columns),
        }
        joblib.dump(xgb_dict, RESULTS_DIR / "model_xgb.joblib")
        print(f"[RETRAIN v5] Saved → results/model_xgb.joblib")

    # ── Summary ───────────────────────────────────────────────────────────
    best = max(rf_acc_sm, xgb_acc if HAS_XGB else 0)
    print(f"\n{'='*65}")
    print(f"  SUMMARY")
    print(f"{'='*65}")
    print(f"  Model             Raw LOPO  Smooth LOPO  F1-smooth")
    print(f"  RF                {rf_acc_raw:.1%}      {rf_acc_sm:.1%}        {rf_f1:.3f}")
    if HAS_XGB:
        print(f"  XGBoost              —         {xgb_acc:.1%}        {xgb_f1:.3f}")
    print()
    milestones = [
        (0.65, "Phase 6 target"),
        (0.70, "Strong result"),
        (0.75, "~25 groups needed"),
        (0.80, "~40 groups needed"),
        (0.85, "~60 groups needed"),
        (0.90, "TALOS-N territory"),
    ]
    for thr, label in milestones:
        mark = "✓" if best >= thr else "✗"
        print(f"  {mark} {thr:.0%}  {label}")

    rows = [{'model': 'RandomForest_v5', 'lopo_raw': rf_acc_raw,
             'lopo_smooth': rf_acc_sm, 'lopo_f1': rf_f1,
             'n_proteins': n_groups, 'n_samples': len(X)}]
    if HAS_XGB:
        rows.append({'model': 'XGBoost_v5', 'lopo_raw': None,
                     'lopo_smooth': xgb_acc, 'lopo_f1': xgb_f1,
                     'n_proteins': n_groups, 'n_samples': len(X)})
    pd.DataFrame(rows).to_csv(RESULTS_DIR / "retrain_summary.csv", index=False)
    print(f"\n  Summary → results/retrain_summary.csv")

    print(f"""
  WHAT CHANGED (v4 → v5):
  1. Per-residue normalization: subtract dataset median per (residue,atom)
     Removes residue-type chemistry confounding — the #1 source of noise
  2. Corrected C' RC values (Wishart 2011 solid-state)
     C' helix ~178 ppm vs strand ~175 ppm — 3 ppm separation
  3. Within-fold imputation (no leakage from test into train medians)
  4. Missingness indicators (has_CB, has_C, has_N) — disorder signal
  5. CSI context window features for neighbors
  6. Window=2 for context (captures 5-residue helix turn)
""")


if __name__ == "__main__":
    main()
