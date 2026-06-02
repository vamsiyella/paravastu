"""
validate_pairs.py — Audit BMRB/PDB pairs for training data quality.

Run this LOCALLY (not in sandbox) before adding entries to SOLID_STATE_ENTRIES:
    conda activate nmr
    cd C:\\Users\\vamsi\\.vscode\\paravastu
    python src/validate_pairs.py

For each candidate pair it checks:
  - BMRB fetch succeeds and is valid NMR-STAR
  - Residue count match between BMRB sequence and PDB first chain
  - CA and CB shift counts (need both for ML)
  - SS annotation presence (_struct_conf + _struct_sheet_range in CIF)
  - Outputs a PASS/WARN/FAIL verdict with reason

Pairs that PASS are safe to add to SOLID_STATE_ENTRIES in pipeline.py.
Pairs that WARN need manual review (slight residue mismatch or limited CB).
Pairs that FAIL should be excluded.
"""

import sys
import requests
import pynmrstar
from pathlib import Path

SRC_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SRC_DIR))

# ---------------------------------------------------------------------------
# All candidate pairs to evaluate
# Already in pipeline: marked with (*)
# ---------------------------------------------------------------------------
CANDIDATES = [
    # (bmrb_id, pdb_id, label)
    # --- Already in pipeline ---
    (17561, "2LBH",  "EETI-II knottin (*)"),
    (16318, "2JWU",  "Ubiquitin microcrystals (*)"),
    (15380, "1LY2",  "GB1 crystal (*)"),
    (16299, "2JSV",  "Thioredoxin (*)"),
    (5969,  "1H4W",  "BPTI (*)"),
    # --- Flagged as problematic in current set ---
    (15409, "2KIB",  "HET-s prion [KNOWN MISMATCH]"),
    (15865, "1M8M",  "SH3 domain [KNOWN MISMATCH]"),
    (6838,  "1YMZ",  "fd coat protein [KNOWN MISMATCH]"),
    (17557, "2KSF",  "M2 proton channel [MEMBRANE - RISKY]"),
    (17948, "2NUZ",  "Calmodulin (*)"),
    # --- New candidates from recommended list ---
    (25123, "1UBQ",  "Ubiquitin MPD crystal"),
    (15156, "2LGI",  "GB1 MAS"),
    (15283, "2OED",  "GB3"),
    (17937, "2M02",  "CAP-Gly domain"),
    (16327, "3L9H",  "DsbA oxidized"),
    (16565, "1TPH",  "TIM apo"),
    (18509, "2C9V",  "Cu,Zn SOD1"),
    (27589, "3NOU",  "Asparaginase II"),
    (18396, "3L9H",  "DsbA nanocrystals"),  # same PDB, different dataset
    (19025, "2M02",  "CAP-Gly 19.9T dataset"),
]

CACHE_DIR = Path(__file__).resolve().parent.parent / "data" / "val_cache"
CACHE_DIR.mkdir(parents=True, exist_ok=True)


def fetch_bmrb(bmrb_id: int) -> str | None:
    for url in [
        f"https://bmrb.io/ftp/pub/bmrb/entry_directories/bmr{bmrb_id}/bmr{bmrb_id}_3.str",
        f"https://bmrb.io/ftp/pub/bmrb/entry_directories/bmr{bmrb_id}/bmr{bmrb_id}.str",
    ]:
        cached = CACHE_DIR / f"bmr{bmrb_id}_3.str"
        if cached.exists():
            text = cached.read_text(errors='replace')
            if text.strip().startswith("data_"):
                return text
            cached.unlink()
        try:
            r = requests.get(url, timeout=25)
            if r.ok and r.text.strip().startswith("data_"):
                cached.write_text(r.text)
                return r.text
        except Exception:
            continue
    return None


def parse_bmrb(text: str) -> tuple[int, int, int]:
    """Returns (n_residues, ca_count, cb_count)."""
    try:
        entry = pynmrstar.Entry.from_string(text)
        loops = entry.get_loops_by_category('Atom_chem_shift')
        seq_ids, ca, cb = set(), 0, 0
        for loop in loops:
            tags = loop.tags
            if 'Seq_ID' not in tags or 'Atom_ID' not in tags:
                continue
            si = tags.index('Seq_ID')
            ai = tags.index('Atom_ID')
            for row in loop.data:
                try:
                    seq_ids.add(int(row[si]))
                    atom = str(row[ai]).strip()
                    if atom == 'CA':
                        ca += 1
                    elif atom == 'CB':
                        cb += 1
                except Exception:
                    pass
        return len(seq_ids), ca, cb
    except Exception:
        return 0, 0, 0


def fetch_cif(pdb_id: str) -> str | None:
    cached = CACHE_DIR / f"{pdb_id}.cif"
    if cached.exists():
        return cached.read_text(errors='replace')
    try:
        r = requests.get(f"https://files.rcsb.org/download/{pdb_id}.cif", timeout=25)
        if r.ok:
            cached.write_text(r.text)
            return r.text
    except Exception:
        pass
    return None


def parse_cif(cif_text: str) -> tuple[int, int, int]:
    """Returns (first_chain_residues, helix_segments, sheet_segments)."""
    # --- Count residues in first chain from _atom_site ---
    residues: set = set()
    atom_cols: list[str] = []
    in_atom = False
    for line in cif_text.splitlines():
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
                    chain = parts[ci['label_asym_id']]
                    seq   = int(parts[ci['label_seq_id']])
                    residues.add((chain, seq))
            except Exception:
                pass

    first_chain = sorted(set(c for c, _ in residues))[0] if residues else None
    n_chain_res = len(set(s for c, s in residues if c == first_chain)) if first_chain else 0

    # --- Count SS annotations ---
    helix_segs = 0
    sheet_segs = 0
    conf_cols: list[str] = []
    sheet_cols: list[str] = []
    in_conf = False
    in_sheet = False
    for line in cif_text.splitlines():
        s = line.strip()
        if s.startswith('_struct_conf.'):
            conf_cols.append(s.split('.')[1].split()[0])
            in_conf = True
            in_sheet = False
            continue
        if s.startswith('_struct_sheet_range.'):
            sheet_cols.append(s.split('.')[1].split()[0])
            in_sheet = True
            in_conf = False
            continue
        if (s.startswith('_') and not s.startswith('_struct_conf') and not s.startswith('_struct_sheet')) or s == 'loop_':
            in_conf = False
            in_sheet = False
        if in_conf and conf_cols and s and not s.startswith('_') and not s.startswith('#'):
            if 'HELX' in s.upper():
                helix_segs += 1
        if in_sheet and sheet_cols and s and not s.startswith('_') and not s.startswith('#'):
            parts = s.split()
            if parts and not parts[0].startswith('_') and parts[0] != 'loop_':
                sheet_segs += 1

    return n_chain_res, helix_segs, sheet_segs


def verdict(bmrb_res, ca, cb, pdb_res, h_segs, s_segs):
    issues = []
    delta = abs(bmrb_res - pdb_res)

    if bmrb_res == 0:
        return "FAIL", "BMRB parse error"
    if pdb_res == 0:
        return "FAIL", "PDB parse error"
    if delta > 15:
        return "FAIL", f"Residue mismatch too large (Δ={delta})"
    if h_segs + s_segs == 0:
        return "FAIL", "No SS annotations in CIF"
    if cb == 0:
        return "FAIL", "No CB shifts deposited — cannot train CA-CB features"
    if ca < 20:
        issues.append(f"low CA count ({ca})")
    if delta > 5:
        issues.append(f"residue Δ={delta}")
    if h_segs == 0:
        issues.append("no helix segments")
    if s_segs == 0:
        issues.append("no sheet segments")

    if issues:
        return "WARN", "; ".join(issues)
    return "PASS", ""


def main():
    print(f"\n{'BMRB':>6} {'PDB':>6} | {'B_res':>5} {'CA':>5} {'CB':>5} | {'P_res':>5} {'H_seg':>5} {'S_seg':>5} | {'Δres':>4} | {'VER':>5}  Reason")
    print("-" * 100)

    for bmrb_id, pdb_id, label in CANDIDATES:
        # BMRB
        text = fetch_bmrb(bmrb_id)
        if text is None:
            print(f"{bmrb_id:>6} {pdb_id:>6} | BMRB FETCH FAILED  — {label}")
            continue
        bmrb_res, ca, cb = parse_bmrb(text)

        # CIF
        cif = fetch_cif(pdb_id)
        if cif is None:
            print(f"{bmrb_id:>6} {pdb_id:>6} | {bmrb_res:>5} {ca:>5} {cb:>5} | PDB FETCH FAILED  — {label}")
            continue
        pdb_res, h_segs, s_segs = parse_cif(cif)

        ver, reason = verdict(bmrb_res, ca, cb, pdb_res, h_segs, s_segs)
        delta = abs(bmrb_res - pdb_res)
        flag = {"PASS": "✓", "WARN": "~", "FAIL": "✗"}[ver]
        print(f"{bmrb_id:>6} {pdb_id:>6} | {bmrb_res:>5} {ca:>5} {cb:>5} | {pdb_res:>5} {h_segs:>5} {s_segs:>5} | {delta:>4} | {flag} {ver:<4}  {reason}  [{label}]")


if __name__ == "__main__":
    main()
