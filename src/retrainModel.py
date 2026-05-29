import pandas as pd
import sys
sys.path.insert(0, 'src')
from pathlib import Path

cache_dir = Path("data/batch_cache")

# Only load files that have a PDB suffix (e.g. bmr17561_2LBH_raw.csv)
# These are the ones produced by run_batch_with_mapping.py with real DSSP labels
# The plain bmr*_raw.csv files are old duplicates without reliable SS labels
good_files = [f for f in cache_dir.glob("bmr*_raw.csv")
              if len(f.stem.split('_')) >= 3]  # has PDB suffix

print(f"Loading {len(good_files)} labeled batch files:")
for f in sorted(good_files):
    print(f"  {f.name}")

all_data = []
for f in sorted(good_files):
    df = pd.read_csv(f)
    source = f.stem  # e.g. bmr17561_2LBH_raw
    df['source'] = source
    all_data.append(df)

merged = pd.concat(all_data, ignore_index=True)
print(f"\nTotal shifts loaded: {len(merged):,}")
print(f"SS distribution: {merged['ss_class'].value_counts().to_dict()}")

# Drop unknown — these can't train the model
labeled = merged[merged['ss_class'] != 'unknown'].copy()
print(f"Labeled shifts (excluding unknown): {len(labeled):,}")
print(f"Labeled SS distribution: {labeled['ss_class'].value_counts().to_dict()}")

# Patch build_feature_matrix to not require all atoms (pivot drops rows with missing atoms)
# We'll do a looser pivot that keeps residues even if only CA is observed
def build_loose_feature_matrix(df, window_size=1):
    """
    Build feature matrix that keeps all labeled residues,
    even if only some atom types are observed (fills missing with NaN).
    """
    import numpy as np

    RANDOM_COIL = {
        'A': {'CA': 52.3, 'CB': 19.1, 'N': 123.8},
        'R': {'CA': 56.1, 'CB': 30.9, 'N': 120.5},
        'N': {'CA': 53.1, 'CB': 38.9, 'N': 118.7},
        'D': {'CA': 54.2, 'CB': 41.1, 'N': 120.4},
        'C': {'CA': 58.2, 'CB': 28.0, 'N': 118.8},
        'Q': {'CA': 55.8, 'CB': 29.4, 'N': 119.8},
        'E': {'CA': 56.6, 'CB': 30.2, 'N': 120.2},
        'G': {'CA': 45.1, 'CB': None, 'N': 108.8},
        'H': {'CA': 55.9, 'CB': 29.4, 'N': 118.2},
        'I': {'CA': 61.1, 'CB': 38.8, 'N': 120.9},
        'L': {'CA': 55.2, 'CB': 42.4, 'N': 121.8},
        'K': {'CA': 56.5, 'CB': 33.1, 'N': 120.4},
        'M': {'CA': 55.6, 'CB': 33.1, 'N': 119.6},
        'F': {'CA': 57.7, 'CB': 39.6, 'N': 120.3},
        'P': {'CA': 63.3, 'CB': 32.1, 'N': 136.5},
        'S': {'CA': 58.3, 'CB': 63.8, 'N': 115.7},
        'T': {'CA': 61.8, 'CB': 69.8, 'N': 113.6},
        'W': {'CA': 57.5, 'CB': 29.9, 'N': 121.3},
        'Y': {'CA': 57.9, 'CB': 38.8, 'N': 120.3},
        'V': {'CA': 62.2, 'CB': 32.9, 'N': 119.9},
    }

    # Pivot: one row per (source, seq_id), columns = atom shifts
    pivot = df.pivot_table(
        index=['source', 'seq_id', 'residue'],
        columns='atom',
        values='shift',
        aggfunc='mean'
    ).reset_index()

    # Get SS label per residue
    ss_per = (
        df[df['ss_class'] != 'unknown']
        .groupby(['source', 'seq_id'])['ss_class']
        .agg(lambda x: x.mode()[0])
        .reset_index()
    )

    pivot = pivot.merge(ss_per, on=['source', 'seq_id'], how='inner')
    pivot = pivot.sort_values(['source', 'seq_id']).reset_index(drop=True)

    feature_cols = []

    # Raw shift features — keep all that exist, NaN for missing
    for atom in ['CA', 'CB', 'N', 'H', 'HA']:
        col = f'shift_{atom}'
        pivot[col] = pivot[atom] if atom in pivot.columns else np.nan
        feature_cols.append(col)

    # RC deviations
    for atom in ['CA', 'CB', 'N']:
        col = f'dev_{atom}'
        shift_col = atom
        if shift_col in pivot.columns:
            pivot[col] = pivot.apply(
                lambda row: (row[shift_col] - RANDOM_COIL.get(row['residue'], {}).get(atom, np.nan))
                if pd.notna(row.get(shift_col)) else np.nan,
                axis=1
            )
        else:
            pivot[col] = np.nan
        feature_cols.append(col)

    # CA-CB difference
    if 'CA' in pivot.columns and 'CB' in pivot.columns:
        pivot['CA_CB_diff'] = pivot['CA'] - pivot['CB']
    else:
        pivot['CA_CB_diff'] = np.nan
    feature_cols.append('CA_CB_diff')

    # Context window — per source to avoid cross-protein contamination
    for offset in [-1, 1]:
        for atom in ['CA', 'CB']:
            base = f'shift_{atom}'
            new_col = f'{atom}_n{offset:+d}'
            pivot[new_col] = (
                pivot.groupby('source')[base]
                .shift(-offset)
            )
            feature_cols.append(new_col)

    X = pivot[feature_cols].copy()
    y = pivot['ss_class']

    valid = y.isin(['helix', 'strand', 'coil'])
    X = X[valid].reset_index(drop=True)
    y = y[valid].reset_index(drop=True)

    print(f"[ML] Feature matrix: {X.shape[0]} samples × {X.shape[1]} features")
    print(f"[ML] Label distribution:\n{y.value_counts()}")
    return X, y

X, y = build_loose_feature_matrix(labeled)

from ml_module import train_random_forest, train_xgboost, evaluate_model, save_model
results_dir = Path("results")

print("\n--- Random Forest ---")
rf_result = train_random_forest(X, y)
save_model(rf_result, results_dir / "model_rf.joblib")
if rf_result.get('feature_importance') is not None:
    print("\nTop 10 features:")
    print(rf_result['feature_importance'].head(10).to_string())

print("\n--- XGBoost ---")
xgb_result = train_xgboost(X, y)
save_model(xgb_result, results_dir / "model_xgb.joblib")

if len(y.unique()) >= 2 and len(X) >= 20:
    print("\n--- Evaluation (RF) ---")
    evaluate_model(rf_result['model'], X, y)

print("\n✓ ML retraining complete")
print(f"  RF accuracy:  {rf_result['cv_scores'].mean():.3f}")
print(f"  XGB accuracy: {xgb_result['cv_scores'].mean():.3f}")