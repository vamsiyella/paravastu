"""
make_reference_db.py - Build reference_db.csv from 8 cached raw CSVs

Usage:
    python make_reference_db.py

Reads all bmr*_raw.csv files from data/batch_cache/
Outputs: data/reference_db.csv
"""

import pandas as pd
from pathlib import Path
import sys

# Find cache directory
cache_dir = Path("data/batch_cache")
if not cache_dir.exists():
    print(f"Error: {cache_dir} not found")
    sys.exit(1)

# Load all raw CSVs
print(f"Loading cached CSVs from {cache_dir}...")
raw_files = list(cache_dir.glob("bmr*_raw.csv"))
print(f"Found {len(raw_files)} files")

all_data = []
for f in raw_files:
    df = pd.read_csv(f)
    all_data.append(df)
    print(f"  {f.name}: {len(df)} rows")

# Combine
merged = pd.concat(all_data, ignore_index=True)
print(f"\nTotal: {len(merged)} shift records")

# Aggregate by (residue, atom, ss_class)
print("Computing statistics...")
stats = (
    merged.groupby(['residue', 'atom', 'ss_class'])['shift']
    .agg(count='count', mean='mean', median='median', std='std', min='min', max='max')
    .reset_index()
)

# Random coil reference
rc = {
    'A': {'CA': 52.3, 'CB': 19.1, 'N': 123.8},
    'C': {'CA': 58.2, 'CB': 28.0, 'N': 118.8},
    'D': {'CA': 54.2, 'CB': 41.1, 'N': 120.4},
    'E': {'CA': 56.6, 'CB': 30.2, 'N': 120.2},
    'F': {'CA': 57.7, 'CB': 39.6, 'N': 120.3},
    'G': {'CA': 45.1, 'CB': None, 'N': 108.8},
    'H': {'CA': 55.9, 'CB': 29.4, 'N': 118.2},
    'I': {'CA': 61.1, 'CB': 38.8, 'N': 120.9},
    'K': {'CA': 56.5, 'CB': 33.1, 'N': 120.4},
    'L': {'CA': 55.2, 'CB': 42.4, 'N': 121.8},
    'M': {'CA': 55.6, 'CB': 33.1, 'N': 119.6},
    'N': {'CA': 53.1, 'CB': 38.9, 'N': 118.7},
    'P': {'CA': 63.3, 'CB': 32.1, 'N': 136.5},
    'Q': {'CA': 55.8, 'CB': 29.4, 'N': 119.8},
    'R': {'CA': 56.1, 'CB': 30.9, 'N': 120.5},
    'S': {'CA': 58.3, 'CB': 63.8, 'N': 115.7},
    'T': {'CA': 61.8, 'CB': 69.8, 'N': 113.6},
    'V': {'CA': 62.2, 'CB': 32.9, 'N': 119.9},
    'W': {'CA': 57.5, 'CB': 29.9, 'N': 121.3},
    'Y': {'CA': 57.9, 'CB': 38.8, 'N': 120.3},
}

# Add RC deviation
stats['rc_deviation'] = stats.apply(
    lambda row: float(row['mean']) - rc.get(row['residue'], {}).get(row['atom'], float('nan'))
    if pd.notna(row['mean']) else None,
    axis=1
)

# Sort
stats = stats.sort_values(['residue', 'atom', 'ss_class']).reset_index(drop=True)

# Write
output = Path("data/reference_db.csv")
output.parent.mkdir(parents=True, exist_ok=True)
stats.to_csv(output, index=False)

print(f"\nDone!")
print(f"  {len(stats)} statistics rows")
print(f"  {output}")
print(f"\nSS breakdown:")
print(stats['ss_class'].value_counts())
