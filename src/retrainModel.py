r"""
retrainModel.py — Retrain ML models on full batch dataset.
Optimized for solid-state NMR data where CB is often absent but HA is available.

Usage:
    conda activate nmr
    cd C:\Users\vamsi\.vscode\paravastu
    python src/retrainModel.py
"""
import pandas as pd
import numpy as np
import sys
from pathlib import Path
sys.path.insert(0, 'src')

# ---------------------------------------------------------------------------
# Load raw shifts — deduplicate by BMRB ID, prefer PDB-named version
# ---------------------------------------------------------------------------
cache_dir = Path("data/batch_cache")

# Group files by BMRB ID
from collections import defaultdict
bmrb_files = defaultdict(list)
for f in sorted(cache_dir.glob("bmr*_raw.csv")):
    bmrb_id = f.stem.split('_')[0]  # "bmr17561"
    bmrb_files[bmrb_id].append(f)

# For each BMRB ID, prefer the file with a PDB suffix (e.g. bmr17561_2LBH_raw.csv)
# If only one file, use it. If multiple, pick the one with the most underscores.
selected_files = []
for bmrb_id, files in sorted(bmrb_files.items()):
    if len(files) == 1:
        selected_files.append(files[0])
    else:
        # Prefer file with PDB suffix (more underscores in stem = has PDB ID)
        best = max(files, key=lambda f: f.stem.count('_'))
        selected_files.append(best)

all_data = []
for f in sorted(selected_files):
    df = pd.read_csv(f)
    labeled_frac = df['ss_class'].isin(['helix', 'strand', 'coil']).mean()
    if labeled_frac >= 0.50:
        df['source'] = f.stem
        all_data.append(df)
        ss = df['ss_class'].value_counts().to_dict()
        print(f"  {f.name}: {len(df)} shifts, {labeled_frac:.0%} labeled | {ss}")
    else:
        print(f"  {f.name}: SKIPPED ({labeled_frac:.0%} labeled)")

if not all_data:
    print("ERROR: No usable cache files found. Run --batch first.")
    sys.exit(1)

merged = pd.concat(all_data, ignore_index=True)
labeled = merged[merged['ss_class'].isin(['helix', 'strand', 'coil'])].copy()

print(f"\nEntries loaded:              {len(all_data)}")
print(f"Total raw shifts:            {len(merged):,}")
print(f"Labeled (helix/strand/coil): {len(labeled):,}")
print(f"SS distribution:             {labeled['ss_class'].value_counts().to_dict()}")

# ---------------------------------------------------------------------------
# Random coil references — extended with HA
# ---------------------------------------------------------------------------
RANDOM_COIL = {
    'A': {'CA': 52.3, 'CB': 19.1, 'N': 123.8, 'H': 8.25, 'HA': 4.32},
    'R': {'CA': 56.1, 'CB': 30.9, 'N': 120.5, 'H': 8.27, 'HA': 4.38},
    'N': {'CA': 53.1, 'CB': 38.9, 'N': 118.7, 'H': 8.38, 'HA': 4.75},
    'D': {'CA': 54.2, 'CB': 41.1, 'N': 120.4, 'H': 8.37, 'HA': 4.76},
    'C': {'CA': 58.2, 'CB': 28.0, 'N': 118.8, 'H': 8.32, 'HA': 4.69},
    'Q': {'CA': 55.8, 'CB': 29.4, 'N': 119.8, 'H': 8.27, 'HA': 4.37},
    'E': {'CA': 56.6, 'CB': 30.2, 'N': 120.2, 'H': 8.36, 'HA': 4.29},
    'G': {'CA': 45.1, 'CB': None,  'N': 108.8, 'H': 8.33, 'HA': 3.97},
    'H': {'CA': 55.9, 'CB': 29.4, 'N': 118.2, 'H': 8.41, 'HA': 4.63},
    'I': {'CA': 61.1, 'CB': 38.8, 'N': 120.9, 'H': 8.22, 'HA': 4.23},
    'L': {'CA': 55.2, 'CB': 42.4, 'N': 121.8, 'H': 8.16, 'HA': 4.38},
    'K': {'CA': 56.5, 'CB': 33.1, 'N': 120.4, 'H': 8.25, 'HA': 4.36},
    'M': {'CA': 55.6, 'CB': 33.1, 'N': 119.6, 'H': 8.28, 'HA': 4.52},
    'F': {'CA': 57.7, 'CB': 39.6, 'N': 120.3, 'H': 8.30, 'HA': 4.66},
    'P': {'CA': 63.3, 'CB': 32.1, 'N': 136.5, 'H': None,  'HA': 4.44},
    'S': {'CA': 58.3, 'CB': 63.8, 'N': 115.7, 'H': 8.31, 'HA': 4.50},
    'T': {'CA': 61.8, 'CB': 69.8, 'N': 113.6, 'H': 8.24, 'HA': 4.35},
    'W': {'CA': 57.5, 'CB': 29.9, 'N': 121.3, 'H': 8.18, 'HA': 4.70},
    'Y': {'CA': 57.9, 'CB': 38.8, 'N': 120.3, 'H': 8.18, 'HA': 4.60},
    'V': {'CA': 62.2, 'CB': 32.9, 'N': 119.9, 'H': 8.19, 'HA': 4.18},
}

# ---------------------------------------------------------------------------
# Feature engineering
# ---------------------------------------------------------------------------
def build_feature_matrix(df, window_size=2):
    pivot = df.pivot_table(
        index=['source', 'seq_id', 'residue'],
        columns='atom',
        values='shift',
        aggfunc='mean'
    ).reset_index()

    ss_per = (
        df[df['ss_class'] != 'unknown']
        .groupby(['source', 'seq_id'])['ss_class']
        .agg(lambda x: x.mode().iloc[0] if len(x) > 0 else 'coil')
        .reset_index()
    )
    pivot = pivot.merge(ss_per, on=['source', 'seq_id'], how='inner')
    pivot = pivot.sort_values(['source', 'seq_id']).reset_index(drop=True)

    feature_cols = []

    for atom in ['CA', 'CB', 'N', 'H', 'HA']:
        col = f'shift_{atom}'
        pivot[col] = pivot[atom] if atom in pivot.columns else np.nan
        feature_cols.append(col)

    for atom in ['CA', 'CB', 'N', 'H', 'HA']:
        col = f'dev_{atom}'
        if atom in pivot.columns:
            pivot[col] = pivot.apply(
                lambda row, a=atom: (
                    row[a] - RANDOM_COIL.get(row['residue'], {}).get(a, np.nan)
                ) if pd.notna(row.get(a)) else np.nan,
                axis=1
            )
        else:
            pivot[col] = np.nan
        feature_cols.append(col)

    if 'CA' in pivot.columns and 'CB' in pivot.columns:
        pivot['CA_CB_diff'] = pivot['CA'] - pivot['CB']
    else:
        pivot['CA_CB_diff'] = np.nan
    feature_cols.append('CA_CB_diff')

    if 'H' in pivot.columns and 'HA' in pivot.columns:
        pivot['H_HA_diff'] = pivot['H'] - pivot['HA']
    else:
        pivot['H_HA_diff'] = np.nan
    feature_cols.append('H_HA_diff')

    for offset in list(range(-window_size, 0)) + list(range(1, window_size + 1)):
        for atom in ['CA', 'N', 'HA']:
            base = f'shift_{atom}'
            new_col = f'{atom}_n{offset:+d}'
            pivot[new_col] = pivot.groupby('source')[base].shift(-offset)
            feature_cols.append(new_col)
        for atom in ['CA', 'HA']:
            dev_base = f'dev_{atom}'
            new_col = f'dev_{atom}_n{offset:+d}'
            pivot[new_col] = pivot.groupby('source')[dev_base].shift(-offset)
            feature_cols.append(new_col)

    X = pivot[feature_cols].copy()
    y = pivot['ss_class']
    valid = y.isin(['helix', 'strand', 'coil'])
    X = X[valid].reset_index(drop=True)
    y = y[valid].reset_index(drop=True)

    print(f"\n[ML] Feature matrix: {X.shape[0]} samples x {X.shape[1]} features")
    print(f"[ML] Label distribution:\n{y.value_counts()}")
    nan_rates = (X.isna().mean() * 100).sort_values(ascending=False).head(5)
    print("[ML] Top NaN rates:", {k: f"{v:.0f}%" for k, v in nan_rates.items()})
    return X, y

# ---------------------------------------------------------------------------
# Train
# ---------------------------------------------------------------------------
X, y = build_feature_matrix(merged, window_size=2)

from ml_module import train_random_forest, train_xgboost, evaluate_model, save_model
results_dir = Path("results")
results_dir.mkdir(exist_ok=True)

print("\n--- Random Forest ---")
rf_result = train_random_forest(X, y, n_estimators=300)
save_model(rf_result, results_dir / "model_rf.joblib")
if rf_result.get('feature_importance') is not None:
    print("\nTop 15 features:")
    print(rf_result['feature_importance'].head(15).to_string())

print("\n--- XGBoost ---")
xgb_result = train_xgboost(X, y)
save_model(xgb_result, results_dir / "model_xgb.joblib")

if len(y.unique()) >= 2 and len(X) >= 20:
    print("\n--- Evaluation (RF) ---")
    evaluate_model(rf_result['model'], X, y)

print("\nML retraining complete")
print(f"  RF  CV accuracy: {rf_result['cv_scores'].mean():.3f} +/- {rf_result['cv_scores'].std():.3f}")
print(f"  XGB CV accuracy: {xgb_result['cv_scores'].mean():.3f} +/- {xgb_result['cv_scores'].std():.3f}")
