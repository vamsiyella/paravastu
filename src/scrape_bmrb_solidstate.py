"""
scrape_bmrb_solidstate.py — Automatically find and validate all solid-state NMR
BMRB entries suitable for ML training.

Run LOCALLY (takes ~10-20 min, makes many network requests):
    conda activate nmr
    cd C:\\Users\\vamsi\\.vscode\\paravastu
    python src/scrape_bmrb_solidstate.py

What it does:
  1. Fetches the BMRB solid-state NMR entry list
  2. For each entry, downloads the .str file and counts CA/CB shifts + residues
  3. Follows the BMRB→PDB cross-reference link
  4. Downloads the PDB CIF and counts chain A residues + SS annotations
  5. Applies all your training criteria and outputs a ranked CSV

Output: data/bmrb_ssNMR_candidates.csv
  - All entries passing strict filters, ranked by quality score
  - Ready to copy into SOLID_STATE_ENTRIES

Strict filters applied:
  - CA shifts >= 20 (enough backbone coverage)
  - CB shifts >= 10 (needed for CA-CB features)
  - Residue count mismatch <= 15
  - SS annotations present in CIF (helix_segs + sheet_segs >= 2)
  - Not already in your current list (skipped but still reported)

Soft filters (reported but not excluded):
  - Helix content (want proteins with helix)
  - Total labeled residues
"""

import sys
import time
import json
import requests
import pynmrstar
import pandas as pd
from pathlib import Path

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
ROOT_DIR  = Path(__file__).resolve().parent.parent
CACHE_DIR = ROOT_DIR / "data" / "scrape_cache"
OUT_CSV   = ROOT_DIR / "data" / "bmrb_ssNMR_candidates.csv"
CACHE_DIR.mkdir(parents=True, exist_ok=True)

# Entries already in your pipeline — still validated but flagged
ALREADY_HAVE = {25123, 15156, 15283, 19025, 15380, 17561, 11512, 16327,
                18024, 18397, 19031, 25005}

# Entries known to be bad — skip entirely
SKIP = {15409, 15865, 6838, 17557, 16299}

REQUEST_DELAY = 0.5   # seconds between requests — be polite to BMRB
MAX_ENTRIES   = 300   # safety cap


# ---------------------------------------------------------------------------
# Network helpers with caching
# ---------------------------------------------------------------------------

def get_json(url: str, cache_key: str) -> dict | None:
    cached = CACHE_DIR / f"{cache_key}.json"
    if cached.exists():
        try:
            return json.loads(cached.read_text())
        except Exception:
            pass
    try:
        time.sleep(REQUEST_DELAY)
        r = requests.get(url, timeout=30)
        if r.ok:
            data = r.json()
            cached.write_text(json.dumps(data))
            return data
    except Exception as e:
        print(f"  [NET] {url}: {e}")
    return None


def get_text(url: str, cache_key: str) -> str | None:
    cached = CACHE_DIR / f"{cache_key}.txt"
    if cached.exists():
        t = cached.read_text(errors='replace')
        if t.strip():
            return t
    try:
        time.sleep(REQUEST_DELAY)
        r = requests.get(url, timeout=45)
        if r.ok and r.text.strip():
            cached.write_text(r.text)
            return r.text
    except Exception as e:
        print(f"  [NET] {url}: {e}")
    return None


# ---------------------------------------------------------------------------
# Step 1: Get BMRB solid-state entry list
# ---------------------------------------------------------------------------

def get_bmrb_solidstate_ids() -> list[int]:
    """
    Fetch the list of solid-state NMR BMRB entry IDs.
    Uses the BMRB API search endpoint.
    """
    print("[STEP 1] Fetching BMRB solid-state entry list...")

    # Method 1: BMRB API search for solid-state entries
    url = "https://api.bmrb.io/v2/search/get_bmrb_ids_from_uniprot/"
    # That's not quite right — use the metadata search instead

    # Method 2: Query BMRB for entries tagged as solid-state
    api_url = "https://api.bmrb.io/v2/list_entries?database=macromolecules"
    data = get_json(api_url, "bmrb_all_entries")

    if data is None:
        # Fallback: hardcoded list from BMRB solid-state page
        print("  [WARN] API failed, using BMRB solid-state page scrape")
        return _get_solidstate_ids_from_page()

    # Filter for solid-state entries using the entry list
    # The API returns a list of entry IDs; we need to check each one
    if isinstance(data, list):
        all_ids = [int(x) for x in data if str(x).isdigit()]
        print(f"  Total BMRB entries: {len(all_ids)}")
        return all_ids[:MAX_ENTRIES]  # limit for testing
    return []


def _get_solidstate_ids_from_page() -> list[int]:
    """Scrape BMRB solid-state NMR page directly."""
    url = "https://bmrb.io/data_library/solidstate.shtml"
    text = get_text(url, "bmrb_solidstate_page")
    if not text:
        # Hard-coded fallback of known solid-state IDs from the BMRB list
        # (curated manually from https://bmrb.io/data_library/solidstate.shtml)
        return [
            5757, 6838, 7111, 7280, 11512, 11741, 11795, 15156, 15283,
            15380, 15409, 15818, 15865, 16299, 16318, 16327, 16391, 16565,
            16566, 17561, 17700, 17937, 17948, 18024, 18396, 18397, 18509,
            18543, 18708, 18731, 19025, 19031, 19754, 19831, 20411, 20706,
            25005, 25123, 25236, 26548, 27402, 27457, 27555, 27556, 27562,
            27589, 27590, 27879, 27959, 28035, 30076, 30088, 50089, 50110,
        ]

    # Parse IDs from the page HTML
    import re
    ids = re.findall(r'bmrbId=(\d+)', text)
    ids += re.findall(r'/summary/\?bmrbId=(\d+)', text)
    ids = list(set(int(x) for x in ids))
    print(f"  Found {len(ids)} IDs from solid-state page")
    return sorted(ids)


# ---------------------------------------------------------------------------
# Step 2: For each entry, fetch BMRB .str and extract metadata
# ---------------------------------------------------------------------------

def parse_bmrb_entry(bmrb_id: int) -> dict:
    """Download and parse BMRB entry. Returns metadata dict."""
    result = {
        'bmrb_id': bmrb_id,
        'n_residues': 0,
        'ca_count': 0,
        'cb_count': 0,
        'n_count': 0,
        'ha_count': 0,
        'pdb_id': None,
        'fetch_ok': False,
    }

    # Try both URL patterns
    text = None
    for url in [
        f"https://bmrb.io/ftp/pub/bmrb/entry_directories/bmr{bmrb_id}/bmr{bmrb_id}_3.str",
        f"https://bmrb.io/ftp/pub/bmrb/entry_directories/bmr{bmrb_id}/bmr{bmrb_id}.str",
    ]:
        text = get_text(url, f"bmr{bmrb_id}")
        if text and text.strip().startswith("data_"):
            break
        text = None

    if not text:
        return result

    result['fetch_ok'] = True

    try:
        entry = pynmrstar.Entry.from_string(text)
    except Exception:
        return result

    # Count shifts by atom type
    seq_ids = set()
    atom_counts = {}
    for loop in entry.get_loops_by_category('Atom_chem_shift'):
        tags = loop.tags
        if 'Seq_ID' not in tags or 'Atom_ID' not in tags:
            continue
        si = tags.index('Seq_ID')
        ai = tags.index('Atom_ID')
        for row in loop.data:
            try:
                seq_ids.add(int(row[si]))
                atom = str(row[ai]).strip().upper()
                atom_counts[atom] = atom_counts.get(atom, 0) + 1
            except Exception:
                pass

    result['n_residues'] = len(seq_ids)
    result['ca_count']   = atom_counts.get('CA', 0)
    result['cb_count']   = atom_counts.get('CB', 0)
    result['n_count']    = atom_counts.get('N', 0)
    result['ha_count']   = atom_counts.get('HA', 0)

    # Try to get linked PDB ID
    pdb_id = _extract_pdb_id(entry)
    result['pdb_id'] = pdb_id

    return result


def _extract_pdb_id(entry) -> str | None:
    """Try to extract PDB ID from BMRB entry metadata."""
    import re
    for sf in entry.frame_list:
        for tag in ['PDB_ID', 'Database_accession_code', 'Related_BMRB_accession_code']:
            try:
                vals = sf.get_tag(tag)
                for v in vals:
                    v = str(v).strip()
                    if re.match(r'^[A-Za-z0-9]{4}$', v) and not v.isdigit():
                        return v.upper()
            except Exception:
                pass
    return None


# ---------------------------------------------------------------------------
# Step 3: Validate PDB match
# ---------------------------------------------------------------------------

def parse_pdb_cif(pdb_id: str) -> dict:
    """Download and parse PDB CIF. Returns SS and residue count info."""
    result = {
        'pdb_id': pdb_id,
        'chain_a_residues': 0,
        'n_chains': 0,
        'helix_segs': 0,
        'sheet_segs': 0,
        'helix_residues': 0,
        'strand_residues': 0,
        'fetch_ok': False,
    }

    url  = f"https://files.rcsb.org/download/{pdb_id}.cif"
    text = get_text(url, f"pdb_{pdb_id}")
    if not text:
        return result

    result['fetch_ok'] = True

    # Count residues per chain
    residues: set = set()
    atom_cols: list = []
    in_atom = False
    for line in text.splitlines():
        s = line.strip()
        if s.startswith('_atom_site.'):
            atom_cols.append(s.split('.')[1].split()[0])
            in_atom = True
            continue
        if in_atom:
            if s.startswith('_') or s.startswith('#') or s == 'loop_':
                in_atom = False
                continue
            if not s:
                continue
            parts = s.split()
            ci = {c: i for i, c in enumerate(atom_cols)}
            try:
                if parts[ci.get('group_PDB', 0)] == 'ATOM':
                    chain  = parts[ci['label_asym_id']]
                    seq_id = int(parts[ci['label_seq_id']])
                    residues.add((chain, seq_id))
            except Exception:
                pass

    chains = sorted(set(c for c, _ in residues))
    result['n_chains'] = len(chains)
    if chains:
        first = chains[0]
        result['chain_a_residues'] = len(set(s for c, s in residues if c == first))

    # Parse SS annotations
    helix_ranges: list[tuple] = []
    sheet_ranges: list[tuple] = []
    conf_cols: list = []
    sheet_cols: list = []
    in_conf = in_sheet = False

    for line in text.splitlines():
        s = line.strip()
        if s.startswith('_struct_conf.'):
            conf_cols.append(s.split('.')[1].split()[0])
            in_conf = True; in_sheet = False; continue
        if s.startswith('_struct_sheet_range.'):
            sheet_cols.append(s.split('.')[1].split()[0])
            in_sheet = True; in_conf = False; continue
        if s.startswith('_') or s == 'loop_':
            in_conf = in_sheet = False

        if in_conf and conf_cols and s and not s.startswith('_') and not s.startswith('#'):
            if 'HELX' in s.upper():
                result['helix_segs'] += 1
                # Try to count helix residues
                parts = s.split()
                ci = {c: i for i, c in enumerate(conf_cols)}
                try:
                    start = int(parts[ci['beg_label_seq_id']])
                    end   = int(parts[ci['end_label_seq_id']])
                    helix_ranges.append((start, end))
                    result['helix_residues'] += (end - start + 1)
                except Exception:
                    pass

        if in_sheet and sheet_cols and s and not s.startswith('_') and not s.startswith('#'):
            parts = s.split()
            if parts and not parts[0].startswith('_') and parts[0] != 'loop_':
                result['sheet_segs'] += 1
                ci = {c: i for i, c in enumerate(sheet_cols)}
                try:
                    start = int(parts[ci['beg_label_seq_id']])
                    end   = int(parts[ci['end_label_seq_id']])
                    result['strand_residues'] += (end - start + 1)
                except Exception:
                    pass

    return result


# ---------------------------------------------------------------------------
# Step 4: Score and filter
# ---------------------------------------------------------------------------

def compute_score(row: dict) -> tuple[bool, str, float]:
    """
    Returns (passes_strict_filter, reason_if_fails, quality_score).
    Higher score = better training data.
    """
    b_res  = row['n_residues']
    ca     = row['ca_count']
    cb     = row['cb_count']
    p_res  = row['chain_a_residues']
    h_segs = row['helix_segs']
    s_segs = row['sheet_segs']
    h_res  = row['helix_residues']
    p_tot  = max(p_res, 1)

    if not row['fetch_ok']:
        return False, "BMRB fetch failed", 0.0
    if not row['pdb_fetch_ok']:
        return False, "PDB fetch failed", 0.0
    if b_res == 0:
        return False, "no residues in BMRB", 0.0
    if p_res == 0:
        return False, "no residues in PDB", 0.0
    if ca < 20:
        return False, f"CA too low ({ca})", 0.0
    if cb < 10:
        return False, f"CB too low ({cb})", 0.0
    if h_segs + s_segs < 2:
        return False, "no SS annotations in CIF", 0.0

    delta = abs(b_res - p_res)
    if delta > 15:
        return False, f"residue mismatch Δ={delta}", 0.0

    # Quality score (higher = better)
    score = 0.0
    score += min(ca, 150) / 150 * 30      # CA coverage (up to 30 pts)
    score += min(cb, 150) / 150 * 30      # CB coverage (up to 30 pts)
    score += (1 - delta / 15) * 10        # residue match quality (up to 10 pts)
    score += min(h_segs, 5) * 3           # helix segments (up to 15 pts)
    score += min(s_segs, 5) * 3           # strand segments (up to 15 pts)
    # Bonus for mixed SS (not all-helix or all-strand)
    h_frac = h_res / p_tot
    if 0.15 < h_frac < 0.85:
        score += 10                        # mixed SS bonus

    return True, "", round(score, 1)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("=" * 65)
    print("  BMRB SOLID-STATE NMR ENTRY SCRAPER")
    print("=" * 65)
    print(f"  Output: {OUT_CSV}")
    print(f"  Cache:  {CACHE_DIR}\n")

    # Step 1: Get entry list
    bmrb_ids = _get_solidstate_ids_from_page()
    print(f"  Candidate entry IDs: {len(bmrb_ids)}")

    # Remove known-bad entries
    bmrb_ids = [x for x in bmrb_ids if x not in SKIP]
    print(f"  After removing known-bad: {len(bmrb_ids)}")

    rows = []
    n_pass = 0
    n_fail = 0
    n_no_pdb = 0

    for i, bmrb_id in enumerate(bmrb_ids):
        print(f"\n[{i+1:3d}/{len(bmrb_ids)}] BMRB {bmrb_id}", end="  ", flush=True)

        # Parse BMRB
        bmrb = parse_bmrb_entry(bmrb_id)
        if not bmrb['fetch_ok']:
            print("BMRB fetch failed")
            n_fail += 1
            continue

        print(f"res={bmrb['n_residues']} CA={bmrb['ca_count']} CB={bmrb['cb_count']}", end="  ")

        # Quick pre-filter before fetching PDB
        if bmrb['ca_count'] < 20 or bmrb['cb_count'] < 10:
            print("→ skip (low CA/CB)")
            n_fail += 1
            continue

        # Get PDB ID
        pdb_id = bmrb['pdb_id']
        if not pdb_id:
            # Try BMRB API for cross-reference
            api = get_json(
                f"https://api.bmrb.io/v2/entry/{bmrb_id}?format=json",
                f"api_{bmrb_id}"
            )
            if api:
                # Search for PDB cross-reference in API response
                import re, json as json_mod
                api_str = json_mod.dumps(api)
                pdb_matches = re.findall(r'"([A-Z0-9]{4})"', api_str)
                # Filter out BMRB-style IDs
                for m in pdb_matches:
                    if not m.isdigit() and len(m) == 4:
                        pdb_id = m
                        break

        if not pdb_id:
            print("→ no PDB link")
            n_no_pdb += 1
            # Still record it — user may want to manually find PDB
            rows.append({
                **bmrb,
                'pdb_id': None,
                'chain_a_residues': 0,
                'n_chains': 0,
                'helix_segs': 0,
                'sheet_segs': 0,
                'helix_residues': 0,
                'strand_residues': 0,
                'pdb_fetch_ok': False,
                'delta_res': 999,
                'helix_pct': 0,
                'strand_pct': 0,
                'passes': False,
                'fail_reason': 'no PDB link',
                'score': 0,
                'already_have': bmrb_id in ALREADY_HAVE,
            })
            continue

        # Parse PDB
        pdb = parse_pdb_cif(pdb_id)
        print(f"PDB {pdb_id} chain_A={pdb['chain_a_residues']} H={pdb['helix_segs']} S={pdb['sheet_segs']}", end="  ")

        delta = abs(bmrb['n_residues'] - pdb['chain_a_residues'])
        h_pct = round(100 * pdb['helix_residues'] / max(pdb['chain_a_residues'], 1), 1)
        s_pct = round(100 * pdb['strand_residues'] / max(pdb['chain_a_residues'], 1), 1)

        row = {
            **bmrb,
            'pdb_id': pdb_id,
            'chain_a_residues': pdb['chain_a_residues'],
            'n_chains': pdb['n_chains'],
            'helix_segs': pdb['helix_segs'],
            'sheet_segs': pdb['sheet_segs'],
            'helix_residues': pdb['helix_residues'],
            'strand_residues': pdb['strand_residues'],
            'pdb_fetch_ok': pdb['fetch_ok'],
            'delta_res': delta,
            'helix_pct': h_pct,
            'strand_pct': s_pct,
        }
        passes, reason, score = compute_score(row)
        row['passes'] = passes
        row['fail_reason'] = reason
        row['score'] = score
        row['already_have'] = bmrb_id in ALREADY_HAVE

        rows.append(row)

        if passes:
            n_pass += 1
            flag = "✓ PASS"
        else:
            n_fail += 1
            flag = f"✗ {reason}"
        print(f"Δ={delta} H={h_pct}% → {flag}")

    # Save results
    df = pd.DataFrame(rows)
    if not df.empty:
        df = df.sort_values(['passes', 'score'], ascending=[False, False]).reset_index(drop=True)
        df.to_csv(OUT_CSV, index=False)

    # Print summary
    print(f"\n\n{'='*65}")
    print(f"  RESULTS SUMMARY")
    print(f"{'='*65}")
    print(f"  Total processed: {len(rows)}")
    print(f"  PASS (strict):   {n_pass}")
    print(f"  FAIL:            {n_fail}")
    print(f"  No PDB link:     {n_no_pdb}")
    print(f"\n  Full results → {OUT_CSV}")

    if not df.empty:
        passing = df[df['passes'] & ~df['already_have']].head(25)
        already = df[df['passes'] & df['already_have']]
        print(f"\n  NEW entries passing all filters ({len(passing)}):")
        print(f"  {'BMRB':>6}  {'PDB':>6}  {'Bres':>4}  {'CA':>4}  {'CB':>4}  {'Δ':>3}  {'H%':>5}  {'S%':>5}  {'Score':>5}")
        print(f"  {'-'*6}  {'-'*6}  {'-'*4}  {'-'*4}  {'-'*4}  {'-'*3}  {'-'*5}  {'-'*5}  {'-'*5}")
        for _, r in passing.iterrows():
            print(f"  {r['bmrb_id']:>6}  {str(r['pdb_id']):>6}  {r['n_residues']:>4}  {r['ca_count']:>4}  "
                  f"{r['cb_count']:>4}  {r['delta_res']:>3}  {r['helix_pct']:>5}  {r['strand_pct']:>5}  {r['score']:>5}")

        print(f"\n  Already-have entries (confirmed valid): {len(already)}")
        for _, r in already.iterrows():
            print(f"    BMRB {r['bmrb_id']}/{r['pdb_id']}: H={r['helix_pct']}% S={r['strand_pct']}% score={r['score']}")

        # Print the SOLID_STATE_ENTRIES snippet
        all_new = df[df['passes'] & ~df['already_have']].head(20)
        print(f"\n\n  SUGGESTED SOLID_STATE_ENTRIES (add to pipeline.py):")
        print(f"  # New entries from automated BMRB scrape")
        for _, r in all_new.iterrows():
            h_note = f"H={r['helix_pct']}%" if r['helix_pct'] > 15 else "strand/coil"
            print(f"  ({int(r['bmrb_id']):>6}, \"{r['pdb_id']}\", \"{r.get('protein_name', 'ssNMR protein')}\"),  "
                  f"# Δ={r['delta_res']} CA={r['ca_count']} CB={r['cb_count']} {h_note}")


if __name__ == "__main__":
    main()
