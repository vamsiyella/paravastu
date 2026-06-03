"""
retrainModel.py — Rebuild ML models with proper protein-aware cross-validation.

Key fixes vs v1:
  1. Leave-One-Protein-Out CV (LOPO-CV) — each fold tests on a held-out protein
     This is the ONLY honest accuracy metric for generalization
  2. No data leakage — train/test split respects protein boundaries  
  3. Deduplication: multiple entries of same protein family are grouped
  4. Reports per-class recall so you can see helix performance specifically
  5. Feature matrix uses window_size=1 with within-protein context only

Run:
    conda activate nmr
    cd C:\\Users\\vamsi\\.vscode\\paravastu
    python src/retrainModel.py
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

# ---------------------------------------------------------------------------
# Protein family grouping — entries of the same protein share a group ID
# This prevents data leakage in LOPO-CV
# ---------------------------------------------------------------------------
PROTEIN_GROUPS = {
    # GB1 family (all nearly identical sequence — treat as ONE group)
    15156: "GB1",
    15283: "GB1",
    15380: "GB1",
    18397: "GB1",
    16873: "GB1",
    30088: "GB1",
    # CAP-Gly family
    19025: "CAP-Gly",
    19031: "CAP-Gly",
    25005: "CAP-Gly",
    17937: "CAP-Gly",
    # Ubiquitin family
    25123: "Ubiquitin",
    11512: "Ubiquitin",
    16318: "Ubiquitin",
    # Unique proteins
    17561: "EETI-II",
    16327: "DsbA",
    18024: "CNBD",
    18543: "DsbA",   # DsbA mutant = same fold
}

# ---------------------------------------------------------------------------
# Random coil references
# ---------------------------------------------------------------------------
RANDOM_COIL = {
    'A': {'CA': 52.3, 'CB': 19.1, 'N': 123.8, 'HA': 4.32},
    'R': {'CA': 56.1, 'CB': 30.9, 'N': 120.5, 'HA': 4.34},
    'N': {'CA': 53.1, 'CB': 38.9, 'N': 118.7, 'HA': 4.75},
    'D': {'CA': 54.2, 'CB': 41.1, 'N': 120.4, 'HA': 4.76},
    'C': {'CA': 58.2, 'CB': 28.0, 'N': 118.8, 'HA': 4.56},
    'Q': {'CA': 55.8, 'CB': 29.4, 'N': 119.8, 'HA': 4.37},
    'E': {'CA': 56.6, 'CB': 30.2, 'N': 120.2, 'HA': 4.29},
    'G': {'CA': 45.1, 'CB': None, 'N': 108.8, 'HA': 3.97},
    'H': {'CA': 55.9, 'CB': 29.4, 'N': 118.2, 'HA': 4.63},
    'I': {'CA': 61.1, 'CB': 38.8, 'N': 120.9, 'HA': 4.17},
    'L': {'CA': 55.2, 'CB': 42.4, 'N': 121.8, 'HA': 4.34},
    'K': {'CA': 56.5, 'CB': 33.1, 'N': 120.4, 'HA': 4.36},
    'M': {'CA': 55.6, 'CB': 33.1, 'N': 119.6, 'HA': 4.52},
    'F': {'CA': 57.7, 'CB': 39.6, 'N': 120.3, 'HA': 4.66},
    'P': {'CA': 63.3, 'CB': 32.1, 'N': 136.5, 'HA': 4.44},
    'S': {'CA': 58.3, 'CB': 63.8, 'N': 115.7, 'HA': 4.47},
    'T': {'CA': 61.8, 'CB': 69.8, 'N': 113.6, 'HA': 4.35},
    'W': {'CA': 57.5, 'CB': 29.9, 'N': 121.3, 'HA': 4.70},
    'Y': {'CA': 57.9, 'CB': 38.8, 'N': 120.3, 'HA': 4.60},
    'V': {'CA': 62.2, 'CB': 32.9, 'N': 119.9, 'HA': 4.18},
}

# ---------------------------------------------------------------------------
# Step 1: Load cached shifts
# ---------------------------------------------------------------------------

def load_all_cached_shifts():
    search_dirs = [CACHE_DIR, RESULTS_DIR]
    all_files: dict[int, Path] = {}

    for d in search_dirs:
        if not d.exists():
            continue
        for f in sorted(d.glob("merged_shifts_*.csv")):
            stem = f.stem
            parts = stem.replace("merged_shifts_", "").split("_")
            bmrb_id = None
            for part in reversed(parts):
                try:
                    bmrb_id = int(part)
                    break
                except ValueError:
                    continue
            if bmrb_id is None:
                continue
            if bmrb_id not in all_files or len(f.stem) > len(all_files[bmrb_id].stem):
                all_files[bmrb_id] = f

    print(f"[RETRAIN] Found {len(all_files)} unique BMRB entries:")

    frames = []
    skipped = []

    for bmrb_id, fpath in sorted(all_files.items()):
        df = pd.read_csv(fpath)
        if 'ss_class' not in df.columns or 'shift' not in df.columns:
            skipped.append((bmrb_id, "missing ss_class or shift column"))
            continue

        labeled = df[df['ss_class'].isin(['helix', 'strand', 'coil'])]
        total_res = df['seq_id'].nunique() if 'seq_id' in df.columns else len(df)
        labeled_res = labeled['seq_id'].nunique() if 'seq_id' in labeled.columns else len(labeled)
        label_frac = labeled_res / max(total_res, 1)

        if label_frac < 0.50:
            skipped.append((bmrb_id, f"only {label_frac:.0%} labeled"))
            continue

        df['source_bmrb'] = bmrb_id
        # Assign protein group for LOPO-CV
        df['protein_group'] = PROTEIN_GROUPS.get(bmrb_id, f"unknown_{bmrb_id}")
        labeled = df[df['ss_class'].isin(['helix', 'strand', 'coil'])]
        frames.append(labeled)

        ss = labeled['ss_class'].value_counts().to_dict()
        grp = PROTEIN_GROUPS.get(bmrb_id, "unique")
        print(f"  BMRB {bmrb_id:>6} [{grp:<12}]: {len(labeled):>5} shifts  "
              f"H={ss.get('helix',0):>4} S={ss.get('strand',0):>4} C={ss.get('coil',0):>4}")

    if skipped:
        print(f"\n  [SKIP] {len(skipped)} entries:")
        for bid, reason in skipped:
            print(f"    BMRB {bid}: {reason}")

    combined = pd.concat(frames, ignore_index=True)

    # Report unique protein groups
    groups = combined['protein_group'].unique()
    print(f"\n[RETRAIN] {len(combined)} shifts from {len(frames)} entries "
          f"({len(groups)} unique protein groups)")
    print(f"  Groups: {sorted(groups)}")
    print(f"  SS: {combined['ss_class'].value_counts().to_dict()}")
    return combined


# ---------------------------------------------------------------------------
# Step 2: Build features
# ---------------------------------------------------------------------------

def build_features(df: pd.DataFrame, window_size: int = 1):
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

    # Raw shifts
    for atom in ['CA', 'CB', 'N', 'H', 'HA', 'C']:
        col = f'shift_{atom}'
        pivot[col] = pivot[atom] if atom in pivot.columns else np.nan
        feature_cols.append(col)

    # RC deviations
    for atom in ['CA', 'CB', 'N', 'HA']:
        col = f'dev_{atom}'
        def get_dev(row, _atom=atom):
            if _atom not in pivot.columns or pd.isna(row.get(_atom, np.nan)):
                return np.nan
            rc = RANDOM_COIL.get(row['residue'], {}).get(_atom)
            return row[_atom] - rc if rc is not None else np.nan
        pivot[col] = pivot.apply(get_dev, axis=1)
        feature_cols.append(col)

    # CA-CB diff
    if 'CA' in pivot.columns and 'CB' in pivot.columns:
        pivot['CA_CB_diff'] = pivot['CA'] - pivot['CB']
        feature_cols.append('CA_CB_diff')

    # Residue one-hot
    for aa in list('ACDEFGHIKLMNPQRSTVWY'):
        col = f'is_{aa}'
        pivot[col] = (pivot['residue'] == aa).astype(int)
        feature_cols.append(col)

    # Context window — within same protein only
    for offset in range(-window_size, window_size + 1):
        if offset == 0:
            continue
        for atom in ['CA', 'CB']:
            base = f'shift_{atom}'
            new_col = f'{atom}_n{offset:+d}'
            pivot[new_col] = pivot.groupby('source_bmrb')[base].shift(-offset)
            feature_cols.append(new_col)

    X = pivot[feature_cols].copy()
    y = pivot['ss_class']
    groups = pivot['protein_group']

    print(f"[RETRAIN] Feature matrix: {X.shape[0]} × {X.shape[1]}")
    print(f"[RETRAIN] Labels: {y.value_counts().to_dict()}")
    return X, y, groups


# ---------------------------------------------------------------------------
# Step 3: Leave-One-Protein-Out CV
# ---------------------------------------------------------------------------

def lopo_cv(model_class, model_kwargs: dict, X: pd.DataFrame, y: pd.Series,
            groups: pd.Series, label: str):
    """
    Leave-One-Protein-OUT cross-validation.
    Each fold: train on all proteins except one, test on that one.
    This is the only honest accuracy metric for generalization.
    """
    lopo = LeaveOneGroupOut()
    unique_groups = sorted(groups.unique())
    n_folds = len(unique_groups)

    print(f"\n[RETRAIN] {label} LOPO-CV ({n_folds} folds = {n_folds} protein groups)")

    all_y_true = []
    all_y_pred = []
    fold_accs = []

    X_filled = X.fillna(X.median(numeric_only=True))

    for fold_idx, (train_idx, test_idx) in enumerate(lopo.split(X_filled, y, groups)):
        test_group = groups.iloc[test_idx[0]]
        X_train = X_filled.iloc[train_idx]
        X_test  = X_filled.iloc[test_idx]
        y_train = y.iloc[train_idx]
        y_test  = y.iloc[test_idx]

        # Class weights on training data only
        classes = sorted(y_train.unique())
        if len(classes) < 2:
            continue
        weights = compute_class_weight('balanced', classes=np.array(classes), y=y_train)
        cw = dict(zip(classes, weights))

        model = model_class(**model_kwargs, class_weight=cw)
        model.fit(X_train, y_train)
        y_pred = model.predict(X_test)

        acc = accuracy_score(y_test, y_pred)
        fold_accs.append(acc)
        all_y_true.extend(y_test)
        all_y_pred.extend(y_pred)

        n_h = (y_test == 'helix').sum()
        n_s = (y_test == 'strand').sum()
        n_c = (y_test == 'coil').sum()
        print(f"  Fold {fold_idx+1:2d} [{test_group:<15}] "
              f"n=H{n_h}/S{n_s}/C{n_c}  acc={acc:.3f}")

    overall_acc = accuracy_score(all_y_true, all_y_pred)
    overall_f1  = f1_score(all_y_true, all_y_pred, average='macro', zero_division=0)
    mean_fold   = np.mean(fold_accs)

    print(f"\n  {label} LOPO results:")
    print(f"  Mean fold accuracy:  {mean_fold:.3f}")
    print(f"  Overall accuracy:    {overall_acc:.3f}")
    print(f"  Macro F1:            {overall_f1:.3f}")
    print(f"\n  Per-class report:")
    print(classification_report(all_y_true, all_y_pred, digits=3, zero_division=0))
    labels = ['coil', 'helix', 'strand']
    cm = confusion_matrix(all_y_true, all_y_pred, labels=labels)
    print(f"  Confusion matrix (rows=true, cols=pred):")
    print(pd.DataFrame(cm, index=labels, columns=labels).to_string())

    return overall_acc, overall_f1, fold_accs


# ---------------------------------------------------------------------------
# Step 4: Final model (train on ALL data)
# ---------------------------------------------------------------------------

def train_final_model(model_class, model_kwargs: dict, X: pd.DataFrame,
                      y: pd.Series, label: str) -> object:
    """Train final model on ALL data with class balancing."""
    X_filled = X.fillna(X.median(numeric_only=True))
    classes = sorted(y.unique())
    weights = compute_class_weight('balanced', classes=np.array(classes), y=y)
    cw = dict(zip(classes, weights))

    model = model_class(**model_kwargs, class_weight=cw)
    model.fit(X_filled, y)
    print(f"[RETRAIN] {label} final model trained on {len(X)} samples")
    return model


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("=" * 65)
    print("  PARAVASTU — ML RETRAINING (v2, LOPO-CV)")
    print("=" * 65)

    df = load_all_cached_shifts()

    n_groups = df['protein_group'].nunique()
    if n_groups < 4:
        print(f"\nWARNING: Only {n_groups} unique protein groups.")
        print("LOPO-CV needs at least 4 groups for meaningful results.")
        print("Add more diverse proteins to SOLID_STATE_ENTRIES and re-run --batch.")

    X, y, groups = build_features(df, window_size=1)
    X_filled_rf = X.fillna(X.median(numeric_only=True))

    print(f"\n{'='*65}")
    print(f"  LEAVE-ONE-PROTEIN-OUT CROSS-VALIDATION")
    print(f"  (This is the honest accuracy — model must generalize to unseen proteins)")
    print(f"{'='*65}")

    # RF LOPO-CV
    rf_acc, rf_f1, rf_folds = lopo_cv(
        RandomForestClassifier,
        dict(n_estimators=300, max_depth=None, min_samples_leaf=2,
             random_state=42, n_jobs=-1),
        X, y, groups, "RandomForest"
    )

    # XGB LOPO-CV
    if HAS_XGB:
        from sklearn.preprocessing import LabelEncoder
        le = LabelEncoder()
        y_enc = pd.Series(le.fit_transform(y), index=y.index)
        groups_enc = groups

        lopo = LeaveOneGroupOut()
        all_true, all_pred = [], []
        fold_accs_xgb = []
        X_filled_xgb = X.fillna(-999)

        print(f"\n[RETRAIN] XGBoost LOPO-CV ({groups.nunique()} folds)")
        for fold_idx, (tr, te) in enumerate(lopo.split(X_filled_xgb, y_enc, groups)):
            test_group = groups.iloc[te[0]]
            X_tr, X_te = X_filled_xgb.iloc[tr], X_filled_xgb.iloc[te]
            y_tr, y_te = y_enc.iloc[tr], y_enc.iloc[te]

            if y_tr.nunique() < 2:
                continue

            counts = np.bincount(y_tr)
            max_c = counts.max()
            sw = np.array([max_c / counts[yi] for yi in y_tr])

            model_xgb = xgb.XGBClassifier(
                n_estimators=300, max_depth=5, learning_rate=0.05,
                subsample=0.8, colsample_bytree=0.8,
                eval_metric='mlogloss', random_state=42, n_jobs=-1, verbosity=0,
            )
            model_xgb.fit(X_tr, y_tr, sample_weight=sw)
            pred = model_xgb.predict(X_te)

            acc = accuracy_score(y_te, pred)
            fold_accs_xgb.append(acc)
            all_true.extend(le.inverse_transform(y_te))
            all_pred.extend(le.inverse_transform(pred))
            print(f"  Fold {fold_idx+1:2d} [{test_group:<15}] acc={acc:.3f}")

        xgb_acc = accuracy_score(all_true, all_pred)
        xgb_f1  = f1_score(all_true, all_pred, average='macro', zero_division=0)
        print(f"\n  XGBoost LOPO results:")
        print(f"  Overall accuracy:  {xgb_acc:.3f}")
        print(f"  Macro F1:          {xgb_f1:.3f}")
        print(classification_report(all_true, all_pred, digits=3, zero_division=0))

    # Train final models on ALL data
    print(f"\n{'='*65}")
    print(f"  FINAL MODELS (trained on all data)")
    print(f"{'='*65}")

    rf_final = train_final_model(
        RandomForestClassifier,
        dict(n_estimators=300, max_depth=None, min_samples_leaf=2,
             random_state=42, n_jobs=-1),
        X, y, "RandomForest"
    )

    fi = pd.Series(rf_final.feature_importances_, index=X.columns).sort_values(ascending=False)
    print("\nTop 15 features (RandomForest):")
    print(fi.head(15).to_string())

    # Save
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    rf_dict = {'model': rf_final, 'feature_importance': fi,
               'lopo_acc': rf_acc, 'lopo_f1': rf_f1}
    joblib.dump(rf_dict, RESULTS_DIR / "model_rf.joblib")
    print(f"\n[RETRAIN] Saved → results/model_rf.joblib")

    if HAS_XGB:
        # Train final XGB
        le2 = LabelEncoder()
        y2 = le2.fit_transform(y)
        X2 = X.fillna(-999)
        counts = np.bincount(y2)
        max_c = counts.max()
        sw2 = np.array([max_c / counts[yi] for yi in y2])
        xgb_final = xgb.XGBClassifier(
            n_estimators=300, max_depth=5, learning_rate=0.05,
            subsample=0.8, colsample_bytree=0.8, eval_metric='mlogloss',
            random_state=42, n_jobs=-1, verbosity=0,
        )
        xgb_final.fit(X2, y2, sample_weight=sw2)
        xgb_dict = {'model': xgb_final, 'label_encoder': le2,
                    'lopo_acc': xgb_acc, 'lopo_f1': xgb_f1}
        joblib.dump(xgb_dict, RESULTS_DIR / "model_xgb.joblib")
        print(f"[RETRAIN] Saved → results/model_xgb.joblib")

    # Summary
    print(f"\n{'='*65}")
    print(f"  SUMMARY")
    print(f"{'='*65}")
    print(f"  RF  LOPO accuracy: {rf_acc:.1%}  F1: {rf_f1:.3f}")
    if HAS_XGB:
        print(f"  XGB LOPO accuracy: {xgb_acc:.1%}  F1: {xgb_f1:.3f}")
    print()

    target = 0.75  # realistic target for LOPO-CV (harder than random split)
    best = max(rf_acc, xgb_acc if HAS_XGB else 0)
    if best >= target:
        print(f"  ✓ TARGET MET ({target:.0%} LOPO accuracy)")
    else:
        print(f"  ✗ {best:.1%} LOPO accuracy. Target: {target:.0%}")
        print(f"  To improve: add more DIVERSE proteins (different folds)")
        print(f"  Current unique folds: {n_groups}")
        print(f"  Need proteins that are NOT: GB1, CAP-Gly, Ubiquitin, DsbA, CNBD, EETI-II")

    # Save summary CSV
    rows = [{'model': 'RandomForest', 'lopo_acc': rf_acc, 'lopo_f1': rf_f1,
             'n_proteins': n_groups, 'n_samples': len(X)}]
    if HAS_XGB:
        rows.append({'model': 'XGBoost', 'lopo_acc': xgb_acc, 'lopo_f1': xgb_f1,
                     'n_proteins': n_groups, 'n_samples': len(X)})
    pd.DataFrame(rows).to_csv(RESULTS_DIR / "retrain_summary.csv", index=False)


if __name__ == "__main__":
    main()