"""
wipe_and_rebuild.py — Nuclear option: delete all cached batch data and retrain from scratch.

Run:
    conda activate nmr
    cd C:\\Users\\vamsi\\.vscode\\paravastu
    python src/wipe_and_rebuild.py

What it deletes:
    data/batch_cache/           <- all per-entry merged_shifts CSVs
    data/reference_db.csv       <- combined stats
    data/batch_log.csv          <- run log
    results/model_rf.joblib     <- trained models
    results/model_xgb.joblib
    results/retrain_summary.csv

What it does NOT touch:
    data/*.cif, data/*.pdb      <- structure files (keep, saves re-downloading)
    data/bmr*.str               <- BMRB NMR-STAR files (keep)
    results/*.png               <- plots
    src/                        <- source code

After running this, re-run the full batch:
    python src/pipeline.py --batch
    python src/retrainModelv5.py
"""

import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
RESULTS = ROOT / "results"
CACHE = DATA / "batch_cache"

DELETE_FILES = [
    DATA / "reference_db.csv",
    DATA / "batch_log.csv",
    RESULTS / "model_rf.joblib",
    RESULTS / "model_xgb.joblib",
    RESULTS / "retrain_summary.csv",
    RESULTS / "combined_stats.csv",
]

print("=" * 60)
print("  WIPE AND REBUILD — Paravastu Pipeline")
print("=" * 60)

# Wipe batch cache directory
if CACHE.exists():
    files = list(CACHE.glob("*"))
    print(f"\n[1] Deleting batch cache ({len(files)} files): {CACHE}")
    for f in files:
        f.unlink()
        print(f"    Deleted: {f.name}")
    print(f"    Done.")
else:
    print(f"\n[1] Batch cache not found (already clean): {CACHE}")

# Wipe individual files
print(f"\n[2] Deleting stale output files:")
for path in DELETE_FILES:
    if path.exists():
        path.unlink()
        print(f"    Deleted: {path}")
    else:
        print(f"    (not found): {path}")

# Also wipe merged_shifts and stats CSVs from results/
print(f"\n[3] Cleaning results/ directory:")
for pattern in ["merged_shifts_*.csv", "stats_*.csv", "predictions_*.csv",
                 "coverage_*.csv", "dssp_*.csv"]:
    for f in RESULTS.glob(pattern):
        f.unlink()
        print(f"    Deleted: {f.name}")

print(f"""
{'='*60}
  DONE. Cache is clean.

  Rebuild steps:
  1. Run batch (takes ~10-20 min):
         python src/pipeline.py --batch

  2. Verify batch log:
         python -c "import pandas as pd; print(pd.read_csv('data/batch_log.csv').to_string())"

  3. Retrain:
         python src/retrainModelv5.py
{'='*60}
""")
