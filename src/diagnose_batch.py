"""
diagnose_batch.py — Run each SOLID_STATE_ENTRIES one at a time and report
what succeeds, what errors, and what hangs.

Run:
    conda activate nmr
    cd C:\\Users\\vamsi\\.vscode\\paravastu
    python src/diagnose_batch.py

This does NOT use --batch. It calls run_pipeline() directly for each entry,
catches all exceptions, and prints a clean summary at the end.
If something hangs, Ctrl+C to break — the summary up to that point still prints.
"""

import sys
import time
import traceback
from pathlib import Path

SRC_DIR  = Path(__file__).resolve().parent
ROOT_DIR = SRC_DIR.parent
sys.path.insert(0, str(SRC_DIR))

from pipeline import run_pipeline, SOLID_STATE_ENTRIES

RESULTS = []

def run_one(bmrb_id, pdb_id, label):
    t0 = time.time()
    try:
        res = run_pipeline(
            bmrb_id=bmrb_id,
            pdb_id=pdb_id,
            run_dssp=True,
            save_results=False,
        )
        elapsed = time.time() - t0

        if not res:
            return {"status": "EMPTY", "elapsed": elapsed, "error": "run_pipeline returned empty dict"}

        merged = res.get('merged_df')
        if merged is None:
            return {"status": "NO_MERGED", "elapsed": elapsed, "error": "no merged_df in results"}

        labeled = merged[merged['ss_class'].isin(['helix','strand','coil'])]
        labeled_frac = len(labeled) / max(len(merged), 1)
        ss_dist = merged['ss_class'].value_counts().to_dict()

        return {
            "status": "OK",
            "elapsed": elapsed,
            "n_shifts": len(merged),
            "labeled_frac": labeled_frac,
            "ss": ss_dist,
            "error": None,
        }

    except KeyboardInterrupt:
        raise
    except Exception as e:
        elapsed = time.time() - t0
        return {
            "status": "ERROR",
            "elapsed": elapsed,
            "error": f"{type(e).__name__}: {e}",
            "traceback": traceback.format_exc(),
        }


def main():
    print(f"\n{'='*65}")
    print(f"  BATCH DIAGNOSTIC — {len(SOLID_STATE_ENTRIES)} entries")
    print(f"{'='*65}")
    print("  Ctrl+C to abort and see partial summary.\n")

    results = []

    for bmrb_id, pdb_id, label in SOLID_STATE_ENTRIES:
        print(f"\n{'─'*65}")
        print(f"  Testing BMRB {bmrb_id} / PDB {pdb_id}  [{label}]")
        print(f"{'─'*65}")

        try:
            result = run_one(bmrb_id, pdb_id, label)
        except KeyboardInterrupt:
            print(f"\n  [ABORTED by user at BMRB {bmrb_id}]")
            results.append({
                "bmrb_id": bmrb_id, "pdb_id": pdb_id, "label": label,
                "status": "ABORTED", "elapsed": 0, "error": "KeyboardInterrupt"
            })
            break

        result["bmrb_id"] = bmrb_id
        result["pdb_id"]  = pdb_id
        result["label"]   = label
        results.append(result)

        if result["status"] == "OK":
            ss = result.get("ss", {})
            print(f"  ✓ OK  {result['elapsed']:.0f}s  |  "
                  f"labeled={result['labeled_frac']:.0%}  |  "
                  f"H={ss.get('helix',0)} S={ss.get('strand',0)} C={ss.get('coil',0)}")
        else:
            print(f"  ✗ {result['status']}  {result['elapsed']:.0f}s")
            print(f"    {result['error']}")
            if result.get('traceback'):
                print("    --- traceback ---")
                for line in result['traceback'].strip().splitlines()[-6:]:
                    print(f"    {line}")

    # ── Summary ───────────────────────────────────────────────────────────
    print(f"\n\n{'='*65}")
    print(f"  SUMMARY")
    print(f"{'='*65}")
    print(f"  {'BMRB':>6}  {'PDB':>6}  {'STATUS':>8}  {'Time':>6}  {'Labeled':>8}  Notes")
    print(f"  {'─'*6}  {'─'*6}  {'─'*8}  {'─'*6}  {'─'*8}  {'─'*30}")

    ok_count   = 0
    fail_count = 0

    for r in results:
        status = r["status"]
        t = f"{r['elapsed']:.0f}s"
        lf = f"{r['labeled_frac']:.0%}" if r.get('labeled_frac') is not None else "—"
        ss = r.get("ss", {})
        notes = ""
        if status == "OK":
            notes = f"H={ss.get('helix',0)} S={ss.get('strand',0)} C={ss.get('coil',0)}"
            ok_count += 1
        else:
            notes = (r.get("error") or "")[:45]
            fail_count += 1
        flag = "✓" if status == "OK" else "✗"
        print(f"  {r['bmrb_id']:>6}  {r['pdb_id']:>6}  {flag} {status:<7}  {t:>6}  {lf:>8}  {notes}")

    print(f"\n  {ok_count} succeeded / {fail_count} failed/aborted")

    fails = [r for r in results if r["status"] != "OK"]
    if fails:
        print(f"\n  FAILING ENTRIES (remove from SOLID_STATE_ENTRIES or fix):")
        for r in fails:
            print(f"    ({r['bmrb_id']}, \"{r['pdb_id']}\", \"{r['label']}\")  ← {r['status']}: {r.get('error','')[:60]}")
    else:
        print(f"\n  All entries passed. Run: python src/pipeline.py --batch")


if __name__ == "__main__":
    main()
