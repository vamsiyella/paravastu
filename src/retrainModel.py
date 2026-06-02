"""
retrainModel.py — Rebuild ML models from the full batch cache.

Run after adding new entries to SOLID_STATE_ENTRIES and running --batch:
    conda activate nmr
    cd C:\\Users\\vamsi\\.vscode\\paravastu
    python src/retrainModel.py

What this does differently from the inline ML in pipeline.py:
  1. Loads ALL per-entry cached shift files from data/batch_cache/
  2. Deduplicates by BMRB ID (preferring PDB-named cache files)
  3. Applies a 50% labeled-residue threshold — drops entries where
     fewer than half the residues got SS labels (bad BMRB/PDB pairing)
  4. Uses SMOTE-style class reweighting to fix strand/helix imbalance
  5. Trains with window_size=1 (not 2 — avoids over-smoothing on small sets)
  6. Adds HA deviation feature (strong helix signal)
  7. Reports per-class precision/recall so you can see strand performance
  8. Saves models to results/ with a training summary CSV

Target: RF and XGBoost both above 90% CV accuracy with recall >= 0.80
        for each of helix, strand, and coil.
"""

import sys
import numpy as np
import pandas as pd
from pathlib import Path

SRC_DIR = Path(__file__).resolve().parent
ROOT_DIR = SRC_DIR.parent
DATA_DIR  = ROOT_DIR / "data"
CACHE_DIR = DATA_DIR / "batch_cache"
RESULTS_DIR = ROOT_DIR / "results"

sys.path.insert(0, str(SRC_DIR))

from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold, cross_val_score, cross_validate
from sklearn.metrics import classification_report, confusion_matrix
from sklearn.utils.class_weight import compute_class_weight
from sklearn.preprocessing import LabelEncoder
import joblib

try:
    import xgboost as xgb
    HAS_XGB = True
except ImportError:
    HAS_XGB = False
    print("[RETRAIN] XGBoost not installed — will train RF only.")


# ---------------------------------------------------------------------------
# Random coil references (Wishart 1995)
# ---------------------------------------------------------------------------
RANDOM_COIL = {
    'A': {'CA': 52.3, 'CB': 19.1, 'N': 123.8, 'H': 8.25, 'HA': 4.32},
    'R': {'CA': 56.1, 'CB': 30.9, 'N': 120.5, 'H': 8.27, 'HA': 4.34},
    'N': {'CA': 53.1, 'CB': 38.9, 'N': 118.7, 'H': 8.38, 'HA': 4.75},
    'D': {'CA': 54.2, 'CB': 41.1, 'N': 120.4, 'H': 8.37, 'HA': 4.76},
    'C': {'CA': 58.2, 'CB': 28.0, 'N': 118.8, 'H': 8.32, 'HA': 4.56},
    'Q': {'CA': 55.8, 'CB': 29.4, 'N': 119.8, 'H': 8.27, 'HA': 4.37},
    'E': {'CA': 56.6, 'CB': 30.2, 'N': 120.2, 'H': 8.36, 'HA': 4.29},
    'G': {'CA': 45.1, 'CB': None, 'N': 108.8, 'H': 8.33, 'HA': 3.97},
    'H': {'CA': 55.9, 'CB': 29.4, 'N': 118.2, 'H': 8.41, 'HA': 4.63},
    'I': {'CA': 61.1, 'CB': 38.8, 'N': 120.9, 'H': 8.22, 'HA': 4.17},
    'L': {'CA': 55.2, 'CB': 42.4, 'N': 121.8, 'H': 8.16, 'HA': 4.34},
    'K': {'CA': 56.5, 'CB': 33.1, 'N': 120.4, 'H': 8.25, 'HA': 4.36},
    'M': {'CA': 55.6, 'CB': 33.1, 'N': 119.6, 'H': 8.28, 'HA': 4.52},
    'F': {'CA': 57.7, 'CB': 39.6, 'N': 120.3, 'H': 8.30, 'HA': 4.66},
    'P': {'CA': 63.3, 'CB': 32.1, 'N': 136.5, 'H': None,  'HA': 4.44},
    'S': {'CA': 58.3, 'CB': 63.8, 'N': 115.7, 'H': 8.31, 'HA': 4.47},
    'T': {'CA': 61.8, 'CB': 69.8, 'N': 113.6, 'H': 8.24, 'HA': 4.35},
    'W': {'CA': 57.5, 'CB': 29.9, 'N': 121.3, 'H': 8.18, 'HA': 4.70},
    'Y': {'CA': 57.9, 'CB': 38.8, 'N': 120.3, 'H': 8.18, 'HA': 4.60},
    'V': {'CA': 62.2, 'CB': 32.9, 'N': 119.9, 'H': 8.19, 'HA': 4.18},
}


# ---------------------------------------------------------------------------
# Step 1: Load and deduplicate cached shift files
# ---------------------------------------------------------------------------

def load_all_cached_shifts() -> pd.DataFrame:
    """
    Load all per-entry merged shift CSV files from data/batch_cache/.

    Deduplication strategy:
      - Cache files are named either  merged_shifts_{bmrb_id}.csv
        or  merged_shifts_{pdb_id}_{bmrb_id}.csv
      - For a given BMRB ID, prefer the PDB-named file if both exist
        (it was produced by a run that had a good BMRB/PDB match).
      - Assign each file a source_bmrb tag so we can track provenance.
    """
    if not CACHE_DIR.exists():
        raise FileNotFoundError(
            f"Batch cache not found: {CACHE_DIR}\n"
            "Run:  python src/pipeline.py --batch\n"
            "first to populate the cache."
        )

    # Also look in results/ for any merged_shifts_*.csv files
    search_dirs = [CACHE_DIR, RESULTS_DIR]
    all_files: dict[int, Path] = {}  # bmrb_id -> best file

    for d in search_dirs:
        for f in sorted(d.glob("merged_shifts_*.csv")):
            # Extract BMRB ID from filename
            stem = f.stem  # e.g. merged_shifts_17561 or merged_shifts_2LBH_17561
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

            # Prefer PDB-named files (have longer stem = more specific)
            if bmrb_id not in all_files or len(f.stem) > len(all_files[bmrb_id].stem):
                all_files[bmrb_id] = f

    if not all_files:
        raise FileNotFoundError(
            "No merged_shifts_*.csv files found in batch_cache/ or results/.\n"
            "Run:  python src/pipeline.py --batch  first."
        )

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

        # 50% labeled threshold — entries with fewer are likely bad pairings
        if label_frac < 0.50:
            skipped.append((bmrb_id, f"only {label_frac:.0%} residues labeled (threshold: 50%)"))
            continue

        df['source_bmrb'] = bmrb_id
        frames.append(labeled)

        ss_dist = labeled['ss_class'].value_counts().to_dict()
        print(f"  BMRB {bmrb_id:>6}: {len(labeled):>5} labeled shifts  "
              f"({label_frac:.0%} labeled)  "
              f"H={ss_dist.get('helix',0)} S={ss_dist.get('strand',0)} C={ss_dist.get('coil',0)}"
              f"  [{fpath.name}]")

    if skipped:
        print(f"\n  [SKIP] {len(skipped)} entries excluded:")
        for bmrb_id, reason in skipped:
            print(f"    BMRB {bmrb_id}: {reason}")

    if not frames:
        raise ValueError("No valid labeled data after filtering. Check your batch cache.")

    combined = pd.concat(frames, ignore_index=True)
    print(f"\n[RETRAIN] Combined dataset: {len(combined)} shift records "
          f"from {combined['source_bmrb'].nunique()} proteins")
    print(f"  SS distribution: {combined['ss_class'].value_counts().to_dict()}")
    return combined


# ---------------------------------------------------------------------------
# Step 2: Build feature matrix
# ---------------------------------------------------------------------------

def build_features(df: pd.DataFrame, window_size: int = 1) -> tuple[pd.DataFrame, pd.Series]:
    """
    Build feature matrix with:
    - Raw shifts: CA, CB, N, HA, H
    - RC deviations: dev_CA, dev_CB, dev_N, dev_HA  (most important features)
    - CA-CB difference (structural fingerprint)
    - Residue one-hot (20 features)
    - Context window shifts ±window_size (CA and CB only)
    """
    # Pivot to one row per residue within each source entry
    # Use (source_bmrb, seq_id) as the residue key
    if 'source_bmrb' not in df.columns:
        df = df.copy()
        df['source_bmrb'] = 0

    pivot = df.pivot_table(
        index=['source_bmrb', 'seq_id', 'residue'],
        columns='atom',
        values='shift',
        aggfunc='mean',
    ).reset_index()

    # Get SS label (most common label per residue)
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
    for atom in ['CA', 'CB', 'N', 'H', 'HA']:
        col = f'shift_{atom}'
        pivot[col] = pivot[atom] if atom in pivot.columns else np.nan
        feature_cols.append(col)

    # RC deviation features (the most informative predictors)
    for atom in ['CA', 'CB', 'N', 'HA']:
        col = f'dev_{atom}'
        def get_dev(row, _atom=atom):
            if _atom not in pivot.columns or pd.isna(row.get(_atom, np.nan)):
                return np.nan
            rc = RANDOM_COIL.get(row['residue'], {}).get(_atom)
            if rc is None:
                return np.nan
            return row[_atom] - rc
        pivot[col] = pivot.apply(get_dev, axis=1)
        feature_cols.append(col)

    # CA - CB difference
    if 'CA' in pivot.columns and 'CB' in pivot.columns:
        pivot['CA_CB_diff'] = pivot['CA'] - pivot['CB']
        feature_cols.append('CA_CB_diff')

    # Residue one-hot
    for aa in list('ACDEFGHIKLMNPQRSTVWY'):
        col = f'is_{aa}'
        pivot[col] = (pivot['residue'] == aa).astype(int)
        feature_cols.append(col)

    # Context window — per-protein sequential shifts
    # Important: only shift WITHIN the same source_bmrb entry
    for offset in range(-window_size, window_size + 1):
        if offset == 0:
            continue
        for atom in ['CA', 'CB']:
            base = f'shift_{atom}'
            new_col = f'{atom}_n{offset:+d}'
            # Group by source_bmrb and shift within group
            pivot[new_col] = (
                pivot.groupby('source_bmrb')[base]
                .shift(-offset)
            )
            feature_cols.append(new_col)

    X = pivot[feature_cols].copy()
    y = pivot['ss_class']

    print(f"[RETRAIN] Feature matrix: {X.shape[0]} samples × {X.shape[1]} features")
    print(f"[RETRAIN] Label distribution:\n{y.value_counts().to_string()}")

    return X, y


# ---------------------------------------------------------------------------
# Step 3: Train with class balancing
# ---------------------------------------------------------------------------

def train_rf(X: pd.DataFrame, y: pd.Series, n_folds: int = 5) -> dict:
    X_filled = X.fillna(X.median(numeric_only=True))

    # Compute class weights to fix strand imbalance
    classes = sorted(y.unique())
    weights = compute_class_weight('balanced', classes=np.array(classes), y=y)
    class_weight_dict = dict(zip(classes, weights))
    print(f"[RETRAIN] RF class weights: {class_weight_dict}")

    model = RandomForestClassifier(
        n_estimators=300,
        max_depth=None,
        min_samples_leaf=2,
        class_weight=class_weight_dict,
        random_state=42,
        n_jobs=-1,
    )

    min_class = y.value_counts().min()
    n_splits = min(n_folds, min_class)
    if n_splits < 2:
        print("[RETRAIN] WARNING: not enough samples per class for CV. Training without CV.")
        model.fit(X_filled, y)
        return {'model': model, 'cv_mean': None, 'cv_std': None,
                'feature_importance': pd.Series(model.feature_importances_, index=X.columns).sort_values(ascending=False)}

    cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=42)
    cv_results = cross_validate(
        model, X_filled, y, cv=cv,
        scoring=['accuracy', 'f1_macro'],
        return_train_score=False,
    )
    cv_acc = cv_results['test_accuracy']
    cv_f1  = cv_results['test_f1_macro']
    print(f"[RETRAIN] RF CV accuracy: {cv_acc.mean():.3f} ± {cv_acc.std():.3f}")
    print(f"[RETRAIN] RF CV F1-macro: {cv_f1.mean():.3f} ± {cv_f1.std():.3f}")

    model.fit(X_filled, y)
    fi = pd.Series(model.feature_importances_, index=X.columns).sort_values(ascending=False)

    return {
        'model': model,
        'cv_mean': float(cv_acc.mean()),
        'cv_std':  float(cv_acc.std()),
        'f1_mean': float(cv_f1.mean()),
        'cv_scores': cv_acc,
        'feature_importance': fi,
        'model_type': 'RandomForest',
    }


def train_xgb(X: pd.DataFrame, y: pd.Series, n_folds: int = 5) -> dict:
    if not HAS_XGB:
        print("[RETRAIN] XGBoost not available. Skipping.")
        return {}

    le = LabelEncoder()
    y_enc = le.fit_transform(y)
    X_filled = X.fillna(-999)

    # XGBoost handles class imbalance via scale_pos_weight per class
    # For multiclass: use sample_weight or scale_pos_weight approach
    counts = np.bincount(y_enc)
    max_count = counts.max()
    # scale_pos_weight is for binary; for multiclass use sample_weight
    sample_weights = np.array([max_count / counts[yi] for yi in y_enc])

    model = xgb.XGBClassifier(
        n_estimators=400,
        max_depth=5,
        learning_rate=0.05,
        subsample=0.8,
        colsample_bytree=0.8,
        eval_metric='mlogloss',
        random_state=42,
        n_jobs=-1,
        verbosity=0,
    )

    min_class = min(np.bincount(y_enc))
    n_splits = min(n_folds, min_class)
    if n_splits < 2:
        model.fit(X_filled, y_enc, sample_weight=sample_weights)
        return {'model': model, 'label_encoder': le, 'cv_mean': None, 'cv_std': None}

    cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=42)
    # CV without sample_weight for fair comparison
    cv_results = cross_validate(
        model, X_filled, y_enc, cv=cv,
        scoring=['accuracy', 'f1_macro'],
    )
    cv_acc = cv_results['test_accuracy']
    cv_f1  = cv_results['test_f1_macro']
    print(f"[RETRAIN] XGB CV accuracy: {cv_acc.mean():.3f} ± {cv_acc.std():.3f}")
    print(f"[RETRAIN] XGB CV F1-macro: {cv_f1.mean():.3f} ± {cv_f1.std():.3f}")

    # Final fit with class balancing via sample_weight
    model.fit(X_filled, y_enc, sample_weight=sample_weights)

    return {
        'model': model,
        'label_encoder': le,
        'cv_mean': float(cv_acc.mean()),
        'cv_std':  float(cv_acc.std()),
        'f1_mean': float(cv_f1.mean()),
        'cv_scores': cv_acc,
        'model_type': 'XGBoost',
    }


# ---------------------------------------------------------------------------
# Step 4: Full evaluation with per-class breakdown
# ---------------------------------------------------------------------------

def full_evaluation(rf_result: dict, xgb_result: dict, X: pd.DataFrame, y: pd.Series):
    print("\n" + "=" * 60)
    print("  HELD-OUT EVALUATION (20% test set)")
    print("=" * 60)

    from sklearn.model_selection import train_test_split

    X_filled_rf = X.fillna(X.median(numeric_only=True))
    X_train, X_test, y_train, y_test = train_test_split(
        X_filled_rf, y, test_size=0.2, stratify=y, random_state=99
    )

    # RF
    rf = rf_result['model']
    rf.fit(X_train, y_train)
    y_pred_rf = rf.predict(X_test)
    print("\n--- Random Forest ---")
    print(classification_report(y_test, y_pred_rf, digits=3))
    print("Confusion matrix (rows=true, cols=pred):")
    labels = ['coil', 'helix', 'strand']
    cm = confusion_matrix(y_test, y_pred_rf, labels=labels)
    print(pd.DataFrame(cm, index=labels, columns=labels).to_string())

    # XGB
    if xgb_result and 'model' in xgb_result:
        le = xgb_result['label_encoder']
        X_filled_xgb = X.fillna(-999)
        _, X_test_x, _, y_test_x = train_test_split(
            X_filled_xgb, y, test_size=0.2, stratify=y, random_state=99
        )
        y_enc_test = le.transform(y_test_x)
        y_pred_xgb = le.inverse_transform(xgb_result['model'].predict(X_test_x))
        print("\n--- XGBoost ---")
        print(classification_report(y_test_x, y_pred_xgb, digits=3))

    # Feature importance
    if rf_result.get('feature_importance') is not None:
        print("\n--- Top 15 features (Random Forest) ---")
        print(rf_result['feature_importance'].head(15).to_string())


# ---------------------------------------------------------------------------
# Step 5: Save
# ---------------------------------------------------------------------------

def save_results(rf_result: dict, xgb_result: dict, summary_rows: list):
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    joblib.dump(rf_result, RESULTS_DIR / "model_rf.joblib")
    print(f"\n[RETRAIN] Saved → results/model_rf.joblib")

    if xgb_result and 'model' in xgb_result:
        joblib.dump(xgb_result, RESULTS_DIR / "model_xgb.joblib")
        print(f"[RETRAIN] Saved → results/model_xgb.joblib")

    if summary_rows:
        pd.DataFrame(summary_rows).to_csv(RESULTS_DIR / "retrain_summary.csv", index=False)
        print(f"[RETRAIN] Saved → results/retrain_summary.csv")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("=" * 60)
    print("  PARAVASTU — ML RETRAINING")
    print("=" * 60)

    # 1. Load data
    df = load_all_cached_shifts()

    # 2. Build features
    X, y = build_features(df, window_size=1)

    n_per_class = y.value_counts()
    if n_per_class.min() < 5:
        print(f"\nWARNING: class '{n_per_class.idxmin()}' has only {n_per_class.min()} samples.")
        print("Run --batch with more entries to improve this class before expecting >90% accuracy.")

    # 3. Train
    print("\n--- Training Random Forest ---")
    rf_result = train_rf(X, y)

    print("\n--- Training XGBoost ---")
    xgb_result = train_xgb(X, y)

    # 4. Evaluate
    full_evaluation(rf_result, xgb_result, X, y)

    # 5. Summary
    summary_rows = []
    for name, res in [("RandomForest", rf_result), ("XGBoost", xgb_result)]:
        if res and res.get('cv_mean') is not None:
            summary_rows.append({
                'model': name,
                'cv_accuracy': res['cv_mean'],
                'cv_std': res['cv_std'],
                'f1_macro': res.get('f1_mean'),
                'n_proteins': df['source_bmrb'].nunique() if 'source_bmrb' in df.columns else '?',
                'n_samples': len(X),
                'n_features': X.shape[1],
            })
            target_met = "✓ TARGET MET" if res['cv_mean'] >= 0.90 else f"✗ need {0.90 - res['cv_mean']:.1%} more"
            print(f"\n{name}: {res['cv_mean']:.1%} CV accuracy  {target_met}")

    # 6. Save
    save_results(rf_result, xgb_result, summary_rows)

    print("\n" + "=" * 60)
    print("  DONE")
    print("=" * 60)


if __name__ == "__main__":
    main()
