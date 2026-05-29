"""Diagnose which batch entries have real SS labels vs all 'unknown'."""
import pandas as pd
from pathlib import Path

cache_dir = Path("data/batch_cache")
print(f"\n{'='*70}")
print("BATCH CACHE DIAGNOSTIC")
print(f"{'='*70}\n")

total_shifts = 0
total_labeled = 0
entry_stats = []

for f in sorted(cache_dir.glob("bmr*_raw.csv")):
    bmrb_id = f.stem.replace("bmr", "").replace("_raw", "")
    df = pd.read_csv(f)
    
    n_total = len(df)
    n_labeled = len(df[df['ss_class'] != 'unknown'])
    pct_labeled = 100 * n_labeled / n_total if n_total > 0 else 0
    
    ss_dist = df['ss_class'].value_counts().to_dict()
    helix = ss_dist.get('helix', 0)
    strand = ss_dist.get('strand', 0)
    coil = ss_dist.get('coil', 0)
    unknown = ss_dist.get('unknown', 0)
    
    total_shifts += n_total
    total_labeled += n_labeled
    
    status = "✓ GOOD" if pct_labeled >= 70 else "✗ POOR"
    entry_stats.append({
        'BMRB': bmrb_id,
        'Total': n_total,
        'Labeled': n_labeled,
        'Pct': f"{pct_labeled:.1f}%",
        'Helix': helix,
        'Strand': strand,
        'Coil': coil,
        'Unknown': unknown,
        'Status': status,
    })

df_stats = pd.DataFrame(entry_stats)
print(df_stats.to_string(index=False))

print(f"\n{'='*70}")
print(f"SUMMARY")
print(f"{'='*70}")
print(f"Total shifts across batch:  {total_shifts:,}")
print(f"Shifts with real SS labels: {total_labeled:,}")
print(f"Percentage labeled:         {100*total_labeled/total_shifts:.1f}%")
print(f"\nTarget: ≥70% labeled (you have {100*total_labeled/total_shifts:.1f}%)")

if 100*total_labeled/total_shifts < 50:
    print("\n⚠️  PROBLEM: Most shifts are 'unknown'. Check:")
    print("   1. Are PDB IDs correctly linked in the batch?")
    print("   2. Did DSSP run for each entry (check results/dssp_*.csv)?")
    print("   3. Did sequence alignment succeed?")
    print("\n→ Re-run batch with pdb_map specified:")
    print("   from pipeline import run_batch")
    print("   run_batch([17561, 15409, 16318, ...],")
    print("             pdb_map={17561: '2LBH', 15409: '2KIB', ...})")
elif 100*total_labeled/total_shifts < 70:
    print("\n⚠️  MARGINAL: Some entries lack SS labels. Fix the worst offenders.")
else:
    print("\n✓ GOOD: Enough labeled data. Retraining should improve if you:")
    print("   1. Add more diverse proteins (different folds)")
    print("   2. Increase window_size parameter in build_feature_matrix()")
    print("   3. Try hyperparameter tuning on the models")
