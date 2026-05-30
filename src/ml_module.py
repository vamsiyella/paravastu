"""
ml_module.py — Machine learning models for secondary structure prediction
               from NMR chemical shifts.

Models:
- RandomForest (fast, interpretable, good baseline)
- XGBoost (usually best performance)
- Feature engineering from shift data

Input:  merged DataFrame (shifts + ss_class labels)
Output: trained model, feature importances, cross-validated accuracy
"""

import numpy as np
import pandas as pd
from pathlib import Path
from typing import Optional, Tuple, Dict

from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.model_selection import (
    cross_val_score, StratifiedKFold, train_test_split
)
from sklearn.preprocessing import LabelEncoder, StandardScaler
from sklearn.metrics import classification_report, confusion_matrix
from sklearn.pipeline import Pipeline
import joblib

try:
    import xgboost as xgb
    HAS_XGB = True
except ImportError:
    HAS_XGB = False
    print("[ML] XGBoost not installed. Using GradientBoosting as fallback.")


# ---------------------------------------------------------------------------
# Feature engineering
# ---------------------------------------------------------------------------

BACKBONE_ATOMS = ['CA', 'CB', 'N', 'H', 'C', 'HA']

# Random coil references (for deviation features)
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

AA_PROPERTIES = {
    # residue -> [hydrophobicity, charge, volume, aromaticity]
    'A': [1.8,  0,  88.6, 0], 'R': [-4.5, 1, 173.4, 0],
    'N': [-3.5, 0, 114.1, 0], 'D': [-3.5, -1, 111.1, 0],
    'C': [2.5,  0, 108.5, 0], 'Q': [-3.5, 0, 143.8, 0],
    'E': [-3.5, -1, 138.4, 0], 'G': [-0.4, 0,  60.1, 0],
    'H': [-3.2, 0.5, 153.2, 1], 'I': [4.5, 0, 166.7, 0],
    'L': [3.8,  0, 166.7, 0], 'K': [-3.9, 1, 168.6, 0],
    'M': [1.9,  0, 162.9, 0], 'F': [2.8,  0, 189.9, 1],
    'P': [-1.6, 0, 112.7, 0], 'S': [-0.8, 0,  89.0, 0],
    'T': [-0.7, 0, 116.1, 0], 'W': [-0.9, 0, 227.8, 1],
    'Y': [-1.3, 0, 193.6, 1], 'V': [4.2,  0, 140.0, 0],
}


def build_feature_matrix(
    merged_df: pd.DataFrame,
    window_size: int = 1,
    include_rc_deviation: bool = True,
    include_aa_properties: bool = True,
) -> Tuple[pd.DataFrame, pd.Series]:
    """
    Build feature matrix and labels for ML training.

    Features per residue:
    - Raw shifts: CA, CB, N (or NaN if missing)
    - RC deviations: CA - RC_CA, CB - RC_CB, N - RC_N
    - AA properties: hydrophobicity, charge, volume, aromaticity
    - Residue type: one-hot encoding
    - Context window shifts (neighbors ± window_size)

    Labels: ss_class (helix / strand / coil)

    Returns: (X DataFrame, y Series)
    """
    # Pivot: one row per residue, columns = atom shifts
    pivot = merged_df.pivot_table(
        index=['seq_id', 'residue'],
        columns='atom',
        values='shift',
        aggfunc='mean'
    ).reset_index()

    # Get SS label (take mode per residue)
    ss_per_residue = (
        merged_df[merged_df['ss_class'] != 'unknown']
        .groupby('seq_id')['ss_class']
        .agg(lambda x: x.mode()[0] if not x.empty else 'coil')
        .reset_index()
    )

    df = pivot.merge(ss_per_residue, on='seq_id', how='inner')
    df = df.sort_values('seq_id').reset_index(drop=True)

    feature_cols = []

    # Raw shift features
    for atom in ['CA', 'CB', 'N', 'H', 'HA', 'C']:
        col = f'shift_{atom}'
        if atom in df.columns:
            df[col] = df[atom]
        else:
            df[col] = np.nan
        feature_cols.append(col)

    # RC deviation features
    if include_rc_deviation:
        for atom in ['CA', 'CB', 'N']:
            col = f'dev_{atom}'
            def get_dev(row, atom=atom):
                res = row['residue']
                shift_col = atom
                if shift_col not in row.index or pd.isna(row.get(shift_col)):
                    return np.nan
                rc = RANDOM_COIL.get(res, {}).get(atom)
                if rc is None:
                    return np.nan
                return row[shift_col] - rc
            df[col] = df.apply(get_dev, axis=1)
            feature_cols.append(col)

    # CA - CB difference (structurally informative)
    if 'CA' in df.columns and 'CB' in df.columns:
        df['CA_CB_diff'] = df['CA'] - df['CB']
        feature_cols.append('CA_CB_diff')

    # AA physicochemical properties
    if include_aa_properties:
        prop_names = ['hydrophobicity', 'charge', 'volume', 'aromaticity']
        for i, prop in enumerate(prop_names):
            df[prop] = df['residue'].apply(
                lambda r: AA_PROPERTIES.get(r, [0, 0, 0, 0])[i]
            )
            feature_cols.append(prop)

    # One-hot residue type
    aas = list('ACDEFGHIKLMNPQRSTVWY')
    for aa in aas:
        col = f'is_{aa}'
        df[col] = (df['residue'] == aa).astype(int)
        feature_cols.append(col)

    # Context window (neighboring residue shifts)
    if window_size > 0:
        for offset in range(-window_size, window_size + 1):
            if offset == 0:
                continue
            for atom in ['CA', 'CB']:
                base_col = f'shift_{atom}'
                if base_col not in df.columns:
                    continue
                new_col = f'{atom}_n{offset:+d}'
                df[new_col] = df[base_col].shift(-offset)
                feature_cols.append(new_col)

    X = df[feature_cols].copy()
    y = df['ss_class']

    # Drop rows where label is unknown/missing
    valid = y.isin(['helix', 'strand', 'coil'])
    X = X[valid].reset_index(drop=True)
    y = y[valid].reset_index(drop=True)

    print(f"[ML] Feature matrix: {X.shape[0]} samples × {X.shape[1]} features")
    print(f"[ML] Label distribution:\n{y.value_counts()}")

    return X, y


# ---------------------------------------------------------------------------
# Model training
# ---------------------------------------------------------------------------

def train_random_forest(
    X: pd.DataFrame,
    y: pd.Series,
    n_estimators: int = 200,
    cv_folds: int = 5,
) -> dict:
    """Train RandomForest with cross-validation."""
    model = Pipeline([
        ('rf', RandomForestClassifier(
            n_estimators=n_estimators,
            max_depth=None,
            min_samples_leaf=2,
            class_weight='balanced',
            random_state=42,
            n_jobs=-1,
        ))
    ])

    X_filled = X.fillna(X.mean(numeric_only=True))

    if len(y.unique()) < 2:
        print("[ML] Not enough class diversity for CV. Training on full data.")
        model.fit(X_filled, y)
        return {'model': model, 'cv_scores': None, 'feature_importance': None}

    cv = StratifiedKFold(n_splits=min(cv_folds, y.value_counts().min()), shuffle=True, random_state=42)
    cv_scores = cross_val_score(model, X_filled, y, cv=cv, scoring='accuracy')
    print(f"[ML] RandomForest CV accuracy: {cv_scores.mean():.3f} ± {cv_scores.std():.3f}")

    model.fit(X_filled, y)
    fi = pd.Series(
        model.named_steps['rf'].feature_importances_,
        index=X.columns
    ).sort_values(ascending=False)

    return {
        'model': model,
        'cv_scores': cv_scores,
        'feature_importance': fi,
        'model_type': 'RandomForest',
    }


def train_xgboost(
    X: pd.DataFrame,
    y: pd.Series,
    cv_folds: int = 5,
) -> dict:
    """Train XGBoost classifier."""
    le = LabelEncoder()
    y_enc = le.fit_transform(y)

    X_filled = X.fillna(-999)  # XGBoost handles this naturally

    if HAS_XGB:
        model = xgb.XGBClassifier(
            n_estimators=300,
            max_depth=5,
            learning_rate=0.05,
            subsample=0.8,
            colsample_bytree=0.8,

            eval_metric='mlogloss',
            random_state=42,
            n_jobs=-1,
        )
    else:
        print("[ML] Using GradientBoosting (XGBoost not available)")
        model = GradientBoostingClassifier(
            n_estimators=200,
            max_depth=4,
            learning_rate=0.05,
            random_state=42,
        )

    if len(np.unique(y_enc)) < 2:
        model.fit(X_filled, y_enc)
        return {'model': model, 'label_encoder': le, 'cv_scores': None}

    cv = StratifiedKFold(n_splits=min(cv_folds, min(np.bincount(y_enc))), shuffle=True, random_state=42)
    cv_scores = cross_val_score(model, X_filled, y_enc, cv=cv, scoring='accuracy')
    print(f"[ML] XGBoost CV accuracy: {cv_scores.mean():.3f} ± {cv_scores.std():.3f}")

    model.fit(X_filled, y_enc)

    return {
        'model': model,
        'label_encoder': le,
        'cv_scores': cv_scores,
        'model_type': 'XGBoost' if HAS_XGB else 'GradientBoosting',
    }


# ---------------------------------------------------------------------------
# Evaluation
# ---------------------------------------------------------------------------

def evaluate_model(model, X: pd.DataFrame, y: pd.Series, label_encoder=None) -> dict:
    """Full evaluation on held-out test set."""
    X_filled = X.fillna(X.mean(numeric_only=True))

    X_train, X_test, y_train, y_test = train_test_split(
        X_filled, y, test_size=0.2, stratify=y, random_state=42
    )

    model.fit(X_train, y_train)
    y_pred = model.predict(X_test)

    print("\n[ML] Test set classification report:")
    print(classification_report(y_test, y_pred))

    cm = confusion_matrix(y_test, y_pred)
    print("[ML] Confusion matrix:")
    print(cm)

    return {
        'y_test': y_test,
        'y_pred': y_pred,
        'confusion_matrix': cm,
    }


def save_model(model_dict: dict, path: Path):
    """Save trained model to disk."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    joblib.dump(model_dict, path)
    print(f"[ML] Model saved → {path}")


def load_model(path: Path) -> dict:
    """Load trained model from disk."""
    return joblib.load(path)


# ---------------------------------------------------------------------------
# Full ML pipeline (convenience function)
# ---------------------------------------------------------------------------

def run_ml_pipeline(
    merged_df: pd.DataFrame,
    results_dir: Path,
    window_size: int = 1,
) -> dict:
    """
    End-to-end ML pipeline:
    1. Build features
    2. Train RF + XGBoost
    3. Evaluate
    4. Save models

    Returns dict with both models and evaluation results.
    """
    results_dir = Path(results_dir)

    print("\n" + "="*60)
    print("  ML PIPELINE — Secondary Structure Prediction")
    print("="*60)

    X, y = build_feature_matrix(merged_df, window_size=window_size)

    if len(X) < 10:
        print("[ML] Not enough labeled samples for training. Need at least 10.")
        return {}

    ml_results = {}

    # Random Forest
    print("\n--- Random Forest ---")
    rf_result = train_random_forest(X, y)
    ml_results['rf'] = rf_result
    save_model(rf_result, results_dir / "model_rf.joblib")

    if rf_result.get('feature_importance') is not None:
        print("\nTop 10 features:")
        print(rf_result['feature_importance'].head(10).to_string())

    # XGBoost
    print("\n--- XGBoost ---")
    xgb_result = train_xgboost(X, y)
    ml_results['xgb'] = xgb_result
    save_model(xgb_result, results_dir / "model_xgb.joblib")

    # Evaluate RF on held-out set
    if len(y.unique()) >= 2 and len(X) >= 20:
        print("\n--- Evaluation (RF) ---")
        eval_result = evaluate_model(rf_result['model'], X, y)
        ml_results['evaluation'] = eval_result

    return ml_results
