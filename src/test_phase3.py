"""
test_phase3.py — Phase 3 smoke test

Run this FIRST to verify batch processing works before running the full entry list.
Tests 3 entries: your known-good entry (17561/2LBH) + 2 new ones.

Usage:
    conda activate nmr
    cd C:\\Users\\vamsi\\.vscode\\paravastu
    python src/test_phase3.py
"""

import sys
import os
from pathlib import Path

# Make sure we can import from src/
sys.path.insert(0, str(Path(__file__).resolve().parent))

from batch_module import run_batch, get_entry_list, _process_single_entry


def test_entry_list():
    """Verify the curated entry list loads correctly."""
    print("=" * 60)
    print("  TEST 1: Entry list")
    print("=" * 60)
    df = get_entry_list()
    print(df.to_string(index=False))
    print(f"\n  Total entries: {len(df)}")
    print(f"  Skipped:       {df['skip'].sum()}")
    print(f"  Will process:  {(~df['skip']).sum()}")
    print("  ✓ Entry list OK\n")
    return df


def test_single_entry_known_good():
    """Test with BMRB 17561 / 2LBH — should always work."""
    print("=" * 60)
    print("  TEST 2: Single entry (known-good: 17561 / 2LBH)")
    print("=" * 60)

    repo_root = Path(__file__).resolve().parent.parent
    data_dir = repo_root / "data"

    try:
        df = _process_single_entry(
            bmrb_id=17561,
            pdb_id="2LBH",
            verbose=True,
            data_dir=data_dir,
        )
        if df is None or df.empty:
            print("  ✗ FAILED: returned empty DataFrame")
            return False

        print(f"\n  Shifts returned: {len(df)}")
        print(f"  SS breakdown: {df['ss_class'].value_counts().to_dict()}")
        print(f"  Columns: {list(df.columns)}")
        print(f"  Sample:\n{df.head(5).to_string(index=False)}")

        # Assertions
        assert 'ss_class' in df.columns, "Missing ss_class column"
        assert 'bmrb_id' in df.columns, "Missing bmrb_id column"
        assert len(df) > 100, f"Expected >100 shifts, got {len(df)}"
        assert set(df['ss_class'].unique()) - {'helix', 'strand', 'coil', 'unknown'} == set(), \
            f"Unexpected ss_class values: {df['ss_class'].unique()}"

        print("  ✓ Single entry OK\n")
        return True

    except Exception as e:
        print(f"  ✗ FAILED: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_mini_batch():
    """Run batch with just 3 entries to verify the merge logic."""
    print("=" * 60)
    print("  TEST 3: Mini batch (3 entries)")
    print("=" * 60)

    # Use known-good + 2 more that are relatively simple structures
    test_pairs = [
        (17561, "2LBH"),   # known-good from Phase 1/2
        (16318, "2JWU"),   # ubiquitin — well-characterized, should parse cleanly
        (15380, "1LY2"),   # GB1 — small, well-annotated
    ]

    repo_root = Path(__file__).resolve().parent.parent
    data_dir = repo_root / "data"

    try:
        stats_df = run_batch(
            bmrb_pdb_pairs=test_pairs,
            output_dir=data_dir,
            verbose=True,
        )

        if stats_df.empty:
            print("  ✗ FAILED: empty stats DataFrame")
            return False

        print(f"\n  Stats rows: {len(stats_df)}")
        print(f"  Columns: {list(stats_df.columns)}")

        if 'ss_class' in stats_df.columns:
            print(f"  SS classes: {stats_df['ss_class'].unique().tolist()}")
            ss_counts = stats_df.groupby('ss_class')['count'].sum()
            print(f"  Observations per SS class:\n{ss_counts}")

        # Check we actually have all 3 SS classes if entries processed correctly
        if 'ss_class' in stats_df.columns:
            ss_classes = set(stats_df['ss_class'].unique())
            if {'helix', 'strand', 'coil'}.issubset(ss_classes):
                print("\n  ✓ All 3 SS classes present in reference DB")
            else:
                print(f"\n  ⚠ Only found: {ss_classes} (some entries may have failed)")

        # Check secondary chemical shift effect is visible
        ala_ca = stats_df[(stats_df['residue'] == 'A') & (stats_df['atom'] == 'CA')]
        if len(ala_ca) > 1:
            print("\n  Secondary chemical shift effect (Ala CA):")
            for _, row in ala_ca.iterrows():
                print(f"    [{row['ss_class']:6s}] mean={row['mean']:.2f} ppm  n={row['count']}")

        print("\n  ✓ Mini batch OK")
        return True

    except Exception as e:
        print(f"  ✗ FAILED: {e}")
        import traceback
        traceback.print_exc()
        return False


def main():
    print("\n" + "=" * 60)
    print("  PHASE 3 SMOKE TEST")
    print("=" * 60 + "\n")

    results = {}

    results['entry_list'] = True  # Just prints, doesn't fail
    test_entry_list()

    results['single_entry'] = test_single_entry_known_good()

    if results['single_entry']:
        results['mini_batch'] = test_mini_batch()
    else:
        print("  Skipping mini batch (single entry test failed)")
        results['mini_batch'] = False

    print("\n" + "=" * 60)
    print("  TEST RESULTS")
    print("=" * 60)
    for test, passed in results.items():
        status = "✓ PASS" if passed else "✗ FAIL"
        print(f"  {status}  {test}")

    all_passed = all(results.values())
    print()
    if all_passed:
        print("  All tests passed! Ready for full batch run:")
        print("  python src/pipeline.py --batch")
    else:
        print("  Some tests failed. Fix errors above before running full batch.")
    print("=" * 60 + "\n")


if __name__ == "__main__":
    main()
