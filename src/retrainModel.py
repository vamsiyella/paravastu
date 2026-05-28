# retrain_ml_on_reference_db.py
import pandas as pd
import sys
sys.path.insert(0, 'src')
from ml_module import run_ml_pipeline

# Load all 8 cached raw CSVs (the merged shift data)
from pathlib import Path
cache_dir = Path("data/batch_cache")
all_data = []
for f in cache_dir.glob("bmr*_raw.csv"):
    df = pd.read_csv(f)
    all_data.append(df)

merged = pd.concat(all_data, ignore_index=True)
print(f"Training on {len(merged)} shift records from {len(all_data)} proteins")
print(f"SS distribution: {merged['ss_class'].value_counts().to_dict()}")

# Train ML models on full dataset
results_dir = Path("results")
ml_results = run_ml_pipeline(merged, results_dir, window_size=1)

print("\n✓ ML retraining complete")
print(f"  RF accuracy:  {ml_results['rf']['cv_scores'].mean():.3f}")
print(f"  XGB accuracy: {ml_results['xgb']['cv_scores'].mean():.3f}")