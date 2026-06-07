"""
retrainModel_v3.py — Two improvements:
1. Drop entries that contribute noise (2KIB 6%, 1YMZ 14%, 1M8M all-coil)
2. Use cascaded binary classifiers: helix-vs-rest THEN strand-vs-coil
   This is the standard approach when one class dominates another in shift space
"""
import pandas as pd
import numpy as np
import sys
sys.path.insert(0, 'src')
from pathlib import Path
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score, StratifiedKFold
from sklearn.pipeline import Pipeline
import joblib

cache_dir = Path("data/batch_cache")

# Only use entries with ≥40% labeled AND that have real structural diversity
# Drop: 15409/2KIB (6%), 6838/1YMZ (14%), 15865/1M8M (all-coil CIF)
GOOD_FILES = {
    'bmr15380_1LY2_raw.csv',   # 99% labeled, strand-rich, has CB
    'bmr16299_2JSV_raw.csv',   # 49% labeled, strand, has CB+N
    'bmr16318_2JWU_raw.csv',   # 54% labeled, strand, has CB+HA
    'bmr17557_2KSF_raw.csv',   # 85% labeled, helix-rich, has CB+HA
    'bmr17561_2LBH_raw.csv',   # 93% labeled, mixed, has H+HA+N
}

print("Loading curated entries only:")
all_data = []
for fname in GOOD_FILES:
    f = cache_dir / fname
    if not f.exists():
        print(f"  MISSING: {fname} — skipping")
        continue
    df = pd.read_csv(f)
    df['source'] = f.stem
    pct = 100*(df['ss_class']!='unknown').mean()
    atoms = sorted(df['atom'].unique())[:6]
    print(f"  {fname}: {len(df)} rows, {pct:.0f}% labeled, atoms: {atoms}")
    all_data.append(df)

merged = pd.concat(all_data, ignore_index=True)
labeled = merged[merged['ss_class'] != 'unknown'].copy()
print(f"\nTotal labeled: {len(labeled)}")
print(f"SS distribution: {labeled['ss_class'].value_counts().to_dict()}")

# ── Feature matrix ─────────────────────────────────────────────────────────
RANDOM_COIL = {
    'A':{'CA':52.3,'CB':19.1,'N':123.8,'H':8.25,'HA':4.32},
    'R':{'CA':56.1,'CB':30.9,'N':120.5,'H':8.27,'HA':4.38},
    'N':{'CA':53.1,'CB':38.9,'N':118.7,'H':8.38,'HA':4.75},
    'D':{'CA':54.2,'CB':41.1,'N':120.4,'H':8.37,'HA':4.76},
    'C':{'CA':58.2,'CB':28.0,'N':118.8,'H':8.32,'HA':4.69},
    'Q':{'CA':55.8,'CB':29.4,'N':119.8,'H':8.27,'HA':4.37},
    'E':{'CA':56.6,'CB':30.2,'N':120.2,'H':8.36,'HA':4.29},
    'G':{'CA':45.1,'CB':None,'N':108.8,'H':8.33,'HA':3.97},
    'H':{'CA':55.9,'CB':29.4,'N':118.2,'H':8.41,'HA':4.63},
    'I':{'CA':61.1,'CB':38.8,'N':120.9,'H':8.22,'HA':4.23},
    'L':{'CA':55.2,'CB':42.4,'N':121.8,'H':8.16,'HA':4.38},
    'K':{'CA':56.5,'CB':33.1,'N':120.4,'H':8.25,'HA':4.36},
    'M':{'CA':55.6,'CB':33.1,'N':119.6,'H':8.28,'HA':4.52},
    'F':{'CA':57.7,'CB':39.6,'N':120.3,'H':8.30,'HA':4.66},
    'P':{'CA':63.3,'CB':32.1,'N':136.5,'H':None,'HA':4.44},
    'S':{'CA':58.3,'CB':63.8,'N':115.7,'H':8.31,'HA':4.50},
    'T':{'CA':61.8,'CB':69.8,'N':113.6,'H':8.24,'HA':4.35},
    'W':{'CA':57.5,'CB':29.9,'N':121.3,'H':8.18,'HA':4.70},
    'Y':{'CA':57.9,'CB':38.8,'N':120.3,'H':8.18,'HA':4.60},
    'V':{'CA':62.2,'CB':32.9,'N':119.9,'H':8.19,'HA':4.18},
}

def build_features(df, window=2):
    pivot = df.pivot_table(
        index=['source','seq_id','residue'],
        columns='atom', values='shift', aggfunc='mean'
    ).reset_index()
    ss_per = (
        df[df['ss_class']!='unknown']
        .groupby(['source','seq_id'])['ss_class']
        .agg(lambda x: x.mode()[0]).reset_index()
    )
    pivot = pivot.merge(ss_per, on=['source','seq_id'], how='inner')
    pivot = pivot.sort_values(['source','seq_id']).reset_index(drop=True)
    
    fcols = []
    for atom in ['CA','CB','N','H','HA']:
        col = f'shift_{atom}'
        pivot[col] = pivot[atom] if atom in pivot.columns else np.nan
        fcols.append(col)
    for atom in ['CA','CB','N','H','HA']:
        col = f'dev_{atom}'
        if atom in pivot.columns:
            pivot[col] = pivot.apply(
                lambda r, a=atom: (r[a] - RANDOM_COIL.get(r['residue'],{}).get(a,np.nan))
                if pd.notna(r.get(a)) else np.nan, axis=1)
        else:
            pivot[col] = np.nan
        fcols.append(col)
    if 'CA' in pivot.columns and 'CB' in pivot.columns:
        pivot['CA_CB_diff'] = pivot['CA'] - pivot['CB']
    else:
        pivot['CA_CB_diff'] = np.nan
    fcols.append('CA_CB_diff')
    if 'H' in pivot.columns and 'HA' in pivot.columns:
        pivot['H_HA_diff'] = pivot['H'] - pivot['HA']
    else:
        pivot['H_HA_diff'] = np.nan
    fcols.append('H_HA_diff')
    # Context window
    for offset in list(range(-window,0))+list(range(1,window+1)):
        for atom in ['CA','CB','N','HA']:
            base = f'shift_{atom}'
            col = f'{atom}_n{offset:+d}'
            pivot[col] = pivot.groupby('source')[base].shift(-offset)
            fcols.append(col)
        for atom in ['CA','CB','HA']:
            dev_base = f'dev_{atom}'
            col = f'dev_{atom}_n{offset:+d}'
            pivot[col] = pivot.groupby('source')[dev_base].shift(-offset)
            fcols.append(col)
    X = pivot[fcols]
    y = pivot['ss_class']
    valid = y.isin(['helix','strand','coil'])
    return X[valid].reset_index(drop=True), y[valid].reset_index(drop=True)

X, y = build_features(labeled, window=2)
print(f"\n[ML] Feature matrix: {X.shape[0]} samples × {X.shape[1]} features")
print(f"[ML] Labels: {y.value_counts().to_dict()}")

X_filled = X.fillna(X.mean(numeric_only=True))

# ── Approach 1: Standard 3-class RF (baseline) ───────────────────────────
print("\n--- Standard 3-class RandomForest ---")
rf3 = RandomForestClassifier(n_estimators=300, class_weight='balanced',
                              min_samples_leaf=2, random_state=42, n_jobs=-1)
cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
scores3 = cross_val_score(rf3, X_filled, y, cv=cv, scoring='balanced_accuracy')
print(f"Balanced accuracy CV: {scores3.mean():.3f} ± {scores3.std():.3f}")
rf3.fit(X_filled, y)

# ── Approach 2: Cascaded binary classifiers ───────────────────────────────
print("\n--- Cascaded Binary Classifiers ---")
print("  Stage 1: helix vs (strand+coil)")
y_stage1 = y.map({'helix':'helix','strand':'other','coil':'other'})
rf_helix = RandomForestClassifier(n_estimators=300, class_weight='balanced',
                                   min_samples_leaf=2, random_state=42, n_jobs=-1)
s1_scores = cross_val_score(rf_helix, X_filled, y_stage1, cv=cv, scoring='balanced_accuracy')
print(f"  Stage 1 balanced accuracy: {s1_scores.mean():.3f} ± {s1_scores.std():.3f}")
rf_helix.fit(X_filled, y_stage1)

print("  Stage 2: strand vs coil (on non-helix residues only)")
non_helix_mask = y != 'helix'
X_nh = X_filled[non_helix_mask].reset_index(drop=True)
y_nh = y[non_helix_mask].reset_index(drop=True)
rf_strand = RandomForestClassifier(n_estimators=300, class_weight='balanced',
                                    min_samples_leaf=2, random_state=42, n_jobs=-1)
cv_nh = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
s2_scores = cross_val_score(rf_strand, X_nh, y_nh, cv=cv_nh, scoring='balanced_accuracy')
print(f"  Stage 2 balanced accuracy: {s2_scores.mean():.3f} ± {s2_scores.std():.3f}")
rf_strand.fit(X_nh, y_nh)

# Simulate cascade on full set
y_pred_cascade = []
for i in range(len(X_filled)):
    row = X_filled.iloc[[i]]
    stage1 = rf_helix.predict(row)[0]
    if stage1 == 'helix':
        y_pred_cascade.append('helix')
    else:
        y_pred_cascade.append(rf_strand.predict(row)[0])

from sklearn.metrics import classification_report, balanced_accuracy_score
print("\nCascade full-set report (optimistic — not CV):")
print(classification_report(y, y_pred_cascade))
print(f"Cascaded balanced accuracy: {balanced_accuracy_score(y, y_pred_cascade):.3f}")

# ── Save models ───────────────────────────────────────────────────────────
results_dir = Path("results")
joblib.dump({'model': rf3, 'type': '3class'}, results_dir / "model_rf.joblib")
joblib.dump({'model': rf_helix, 'stage': 1, 
             'model2': rf_strand, 'stage2': 2,
             'type': 'cascade'}, results_dir / "model_cascade.joblib")
print(f"\n✓ Saved model_rf.joblib (3-class) and model_cascade.joblib (cascade)")
print(f"\nSUMMARY:")
print(f"  3-class RF balanced accuracy:  {scores3.mean():.3f}")
print(f"  Cascade Stage 1 (helix):       {s1_scores.mean():.3f}")
print(f"  Cascade Stage 2 (strand/coil): {s2_scores.mean():.3f}")
