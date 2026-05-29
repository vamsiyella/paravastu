"""
retrainModel.py — Retrain ML models on full batch dataset.
Optimized for solid-state NMR data where CB is often absent but HA is available.
Key insight: HA deviation from random coil is highly discriminating for helix vs strand.
"""
import pandas as pd
import numpy as np
import sys
sys.path.insert(0, 'src')
from pathlib import Path

cache_dir = Path("data/batch_cache")

# Only load PDB-suffixed files (real DSSP labels)
good_files = [f for f in cache_dir.glob("bmr*_raw.csv")
              if len(f.stem.split('_')) >= 3]

print(f"Loading {len(good_files)} labeled batch files:")
for f in sorted(good_files):
    df_check = pd.read_csv(f)
    labeled_pct = 100 * (df_check['ss_class'] != 'unknown').mean()
    atoms = sorted(df_check['atom'].unique())[:8]
    print(f"  {f.name}: {len(df_check)} shifts, {labeled_pct:.0f}% labeled, atoms: {atoms}")

all_data = []
for f in sorted(good_files):
    df = pd.read_csv(f)
    df['source'] = f.stem
    all_data.append(df)

merged = pd.concat(all_data, ignore_index=True)
labeled = merged[merged['ss_class'] != 'unknown'].copy()
print(f"\nTotal labeled shifts: {len(labeled):,}")
print(f"SS distribution: {labeled['ss_class'].value_counts().to_dict()}")

# Random coil references — extended with HA
RANDOM_COIL = {
    'A': {'CA': 52.3, 'CB': 19.1, 'N': 123.8, 'H': 8.25, 'HA': 4.32},
    'R': {'CA': 56.1, 'CB': 30.9, 'N': 120.5, 'H': 8.27, 'HA': 4.38},
    'N': {'CA': 53.1, 'CB': 38.9, 'N': 118.7, 'H': 8.38, 'HA': 4.75},
    'D': {'CA': 54.2, 'CB': 41.1, 'N': 120.4, 'H': 8.37, 'HA': 4.76},
    'C': {'CA': 58.2, 'CB': 28.0, 'N': 118.8, 'H': 8.32, 'HA': 4.69},
    'Q': {'CA': 55.8, 'CB': 29.4, 'N': 119.8, 'H': 8.27, 'HA': 4.37},
    'E': {'CA': 56.6, 'CB': 30.2, 'N': 120.2, 'H': 8.36, 'HA': 4.29},
    'G': {'CA': 45.1, 'CB': None, 'N': 108.8, 'H': 8.33, 'HA': 3.97},
    'H': {'CA': 55.9, 'CB': 29.4, 'N': 118.2, 'H': 8.41, 'HA': 4.63},
    'I': {'CA': 61.1, 'CB': 38.8, 'N': 120.9, 'H': 8.22, 'HA': 4.23},
    'L': {'CA': 55.2, 'CB': 42.4, 'N': 121.8, 'H': 8.16, 'HA': 4.38},
    'K': {'CA': 56.5, 'CB': 33.1, 'N': 120.4, 'H': 8.25, 'HA': 4.36},
    'M': {'CA': 55.6, 'CB': 33.1, 'N': 119.6, 'H': 8.28, 'HA': 4.52},
    'F': {'CA': 57.7, 'CB': 39.6, 'N': 120.3, 'H': 8.30, 'HA': 4.66},
    'P': {'CA': 63.3, 'CB': 32.1, 'N': 136.5, 'H': None, 'HA': 4.44},
    'S': {'CA': 58.3, 'CB': 63.8, 'N': 115.7, 'H': 8.31, 'HA': 4.50},
    'T': {'CA': 61.8, 'CB': 69.8, 'N': 113.6, 'H': 8.24, 'HA': 4.35},
    'W': {'CA': 57.5, 'CB': 29.9, 'N': 121.3, 'H': 8.18, 'HA': 4.70},
    'Y': {'CA': 57.9, 'CB': 38.8, 'N': 120.3, 'H': 8.18, 'HA': 4.60},
    'V': {'CA': 62.2, 'CB': 32.9, 'N': 119.9, 'H': 8.19, 'HA': 4.18},
}

def build_feature_matrix(df, window_size=2):
    """
    Build feature matrix optimized for ssNMR data (CA, N, H, HA available; CB often absent).
    Uses HA deviation — one of the best helix/strand discriminators available.
    window_size=2 captures longer-range structural context.
    """
    # Pivot to one row per residue
    pivot = df.pivot_table(
        index=['source', 'seq_id', 'residue'],
        columns='atom',
        values='shift',
        aggfunc='mean'
    ).reset_index()

    # SS label per residue
    ss_per = (
        df[df['ss_class'] != 'unknown']
        .groupby(['source', 'seq_id'])['ss_class']
        .agg(lambda x: x.mode()[0])
        .reset_index()
    )
    pivot = pivot.merge(ss_per, on=['source', 'seq_id'], how='inner')
    pivot = pivot.sort_values(['source', 'seq_id']).reset_index(drop=True)

    feature_cols = []

    # Raw shifts
    for atom in ['CA', 'CB', 'N', 'H', 'HA']:
        col = f'shift_{atom}'
        pivot[col] = pivot[atom] if atom in pivot.columns else np.nan
        feature_cols.append(col)

    # RC deviations — now includes HA which is key for strand detection
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

    # CA-CB difference (when available)
    ca_col = 'CA' if 'CA' in pivot.columns else None
    cb_col = 'CB' if 'CB' in pivot.columns else None
    if ca_col and cb_col:
        pivot['CA_CB_diff'] = pivot[ca_col] - pivot[cb_col]
    else:
        pivot['CA_CB_diff'] = np.nan
    feature_cols.append('CA_CB_diff')

    # H-HA difference (new: captures local environment)
    if 'H' in pivot.columns and 'HA' in pivot.columns:
        pivot['H_HA_diff'] = pivot['H'] - pivot['HA']
    else:
        pivot['H_HA_diff'] = np.nan
    feature_cols.append('H_HA_diff')

    # Context window — group within source to avoid cross-protein bleed
    for offset in list(range(-window_size, 0)) + list(range(1, window_size + 1)):
        for atom in ['CA', 'N', 'HA']:
            base = f'shift_{atom}'
            new_col = f'{atom}_n{offset:+d}'
            pivot[new_col] = pivot.groupby('source')[base].shift(-offset)
            feature_cols.append(new_col)

        # Also include RC deviations in context window for CA and HA
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

    print(f"\n[ML] Feature matrix: {X.shape[0]} samples × {X.shape[1]} features")
    print(f"[ML] Label distribution:\n{y.value_counts()}")
    print(f"[ML] Feature NaN rates (top 10 worst):")
    nan_rates = (X.isna().mean() * 100).sort_values(ascending=False).head(10)
    for feat, rate in nan_rates.items():
        print(f"     {feat}: {rate:.1f}%")
    return X, y

X, y = build_feature_matrix(labeled, window_size=2)

from ml_module import train_random_forest, train_xgboost, evaluate_model, save_model
results_dir = Path("results")

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

print("\n✓ ML retraining complete")
print(f"  RF accuracy:  {rf_result['cv_scores'].mean():.3f} ± {rf_result['cv_scores'].std():.3f}")
print(f"  XGB accuracy: {xgb_result['cv_scores'].mean():.3f} ± {xgb_result['cv_scores'].std():.3f}")