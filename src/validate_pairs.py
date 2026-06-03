"""
validate_pairs.py — Audit BMRB/PDB pairs for training data quality.

Run LOCALLY:
    conda activate nmr
    cd C:\\Users\\vamsi\\.vscode\\paravastu
    python src/validate_pairs.py
"""

import sys
import requests
import pynmrstar
from pathlib import Path

SRC_DIR   = Path(__file__).resolve().parent
ROOT_DIR  = SRC_DIR.parent
CACHE_DIR = ROOT_DIR / "data" / "val_cache"
CACHE_DIR.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# New candidates from ChatGPT list (pre-filtered: removed membrane + all-helix)
# Already confirmed working: 25123/1UBQ, 15156/2LGI, 15283/2OED,
#                            19025/2M02, 15380/1PGB, 15380/2QMT
# Removed: 16299/1XOB (wrong protein — sequences don't match)
# ---------------------------------------------------------------------------
CANDIDATES = [
    # New entries to validate
    (11512, "3ONS",  "Ubiquitin — different dataset, Δ=4"),
    (16327, "1FVK",  "DsbA oxidized — Δ=0, mixed fold"),
    (16565, "7TIM",  "TIM barrel apo — Δ=0, large α/β"),
    (16566, "7TIM",  "TIM barrel +G3P — same PDB, independent shifts"),
    (17700, "1F6M",  "Thioredoxin wt — Δ=3, correct PDB this time"),
    (18024, "2K0G",  "CNBD domain — Δ=0, α/β"),
    (18396, "1A23",  "DsbA reduced nanocrystal — Δ=0"),
    (18397, "1GB1",  "GB1 proton-detected — Δ=0"),
    (18509, "2LU5",  "SOD1 — Δ=0, β-barrel"),
    (18543, "1FVK",  "DsbA C33S mutant — same PDB as 16327"),
    (18708, "2LU5",  "Co-SOD1 — same PDB as 18509, independent shifts"),
    (19031, "2MPX",  "CAP-Gly+EB1 complex — Δ=0"),
    (19831, "2YLK",  "CBM3b cellulose-binding — Δ=0, β-rich"),
    (25005, "2MPX",  "CAP-Gly on microtubule — same PDB as 19031"),
    (27589, "3ECA",  "Asparaginase II — Δ=0, large α/β enzyme"),
    (27590, "3ECA",  "Asparaginase PEG — same PDB, independent shifts"),
]

# ---------------------------------------------------------------------------
# Helpers (same as before)
# ---------------------------------------------------------------------------

def fetch_bmrb(bmrb_id: int) -> str | None:
    cached = CACHE_DIR / f"bmr{bmrb_id}_3.str"
    if cached.exists():
        t = cached.read_text(errors='replace')
        if t.strip().startswith("data_"):
            return t
        cached.unlink()
    for url in [
        f"https://bmrb.io/ftp/pub/bmrb/entry_directories/bmr{bmrb_id}/bmr{bmrb_id}_3.str",
        f"https://bmrb.io/ftp/pub/bmrb/entry_directories/bmr{bmrb_id}/bmr{bmrb_id}.str",
    ]:
        try:
            r = requests.get(url, timeout=25)
            if r.ok and r.text.strip().startswith("data_"):
                cached.write_text(r.text)
                return r.text
        except:
            pass
    return None


def parse_bmrb(text: str) -> tuple[int, int, int]:
    try:
        entry = pynmrstar.Entry.from_string(text)
        loops = entry.get_loops_by_category('Atom_chem_shift')
        seq_ids, ca, cb = set(), 0, 0
        for loop in loops:
            tags = loop.tags
            if 'Seq_ID' not in tags or 'Atom_ID' not in tags:
                continue
            si, ai = tags.index('Seq_ID'), tags.index('Atom_ID')
            for row in loop.data:
                try:
                    seq_ids.add(int(row[si]))
                    a = str(row[ai]).strip()
                    if a == 'CA': ca += 1
                    elif a == 'CB': cb += 1
                except:
                    pass
        return len(seq_ids), ca, cb
    except:
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
    except:
        pass
    return None


def parse_cif(cif_text: str) -> tuple[int, int, int, int]:
    residues: set = set()
    atom_cols: list = []
    in_atom = False
    for line in cif_text.splitlines():
        s = line.strip()
        if s.startswith('_atom_site.'):
            atom_cols.append(s.split('.')[1].split()[0])
            in_atom = True
            continue
        if in_atom:
            if s.startswith('_') or s.startswith('#') or s == 'loop_':
                in_atom = False; continue
            if not s: continue
            parts = s.split()
            ci = {c: i for i, c in enumerate(atom_cols)}
            try:
                if parts[ci.get('group_PDB', 0)] == 'ATOM':
                    residues.add((parts[ci['label_asym_id']], int(parts[ci['label_seq_id']])))
            except:
                pass

    chains = sorted(set(c for c, _ in residues))
    first_chain = chains[0] if chains else None
    n_res = len(set(s for c, s in residues if c == first_chain)) if first_chain else 0
    n_chains = len(chains)

    helix_segs = sheet_segs = 0
    conf_cols: list = []
    sheet_cols: list = []
    in_conf = in_sheet = False
    for line in cif_text.splitlines():
        s = line.strip()
        if s.startswith('_struct_conf.'):
            conf_cols.append(s.split('.')[1].split()[0]); in_conf = True; in_sheet = False; continue
        if s.startswith('_struct_sheet_range.'):
            sheet_cols.append(s.split('.')[1].split()[0]); in_sheet = True; in_conf = False; continue
        if s.startswith('_') or s == 'loop_': in_conf = in_sheet = False
        if in_conf and conf_cols and s and not s.startswith('_') and not s.startswith('#'):
            if 'HELX' in s.upper(): helix_segs += 1
        if in_sheet and sheet_cols and s and not s.startswith('_') and not s.startswith('#'):
            p = s.split()
            if p and not p[0].startswith('_') and p[0] != 'loop_': sheet_segs += 1

    return n_res, n_chains, helix_segs, sheet_segs


def verdict(bmrb_res, ca, cb, pdb_res, n_chains, h_segs, s_segs):
    if bmrb_res == 0: return "FAIL", "BMRB parse error"
    if pdb_res == 0:  return "FAIL", "PDB parse error"
    delta = abs(bmrb_res - pdb_res)
    if h_segs + s_segs == 0: return "FAIL", "No SS annotations in CIF"
    if delta > 15: return "FAIL", f"Residue mismatch (Δ={delta})"
    issues = []
    if cb == 0:     issues.append("CB=0")
    if delta > 5:   issues.append(f"Δ={delta}")
    if ca < 20:     issues.append(f"low CA={ca}")
    if h_segs == 0: issues.append("no helix in CIF")
    if s_segs == 0: issues.append("no strand in CIF")
    return ("WARN" if issues else "PASS"), "; ".join(issues)


def main():
    bmrb_cache: dict[int, tuple] = {}
    prev_bmrb = None

    print(f"\n{'BMRB':>6} {'PDB':>6} | {'B_res':>5} {'CA':>5} {'CB':>5} | {'P_res':>5} {'Nch':>3} {'H':>4} {'S':>4} | {'Δ':>4} | {'VER':>5}  Notes")
    print("-" * 112)

    for bmrb_id, pdb_id, notes in CANDIDATES:
        if bmrb_id != prev_bmrb:
            print()
            prev_bmrb = bmrb_id

        if bmrb_id not in bmrb_cache:
            text = fetch_bmrb(bmrb_id)
            bmrb_cache[bmrb_id] = parse_bmrb(text) if text else (0, 0, 0)
        bmrb_res, ca, cb = bmrb_cache[bmrb_id]

        cif = fetch_cif(pdb_id)
        pdb_res, n_chains, h, s = parse_cif(cif) if cif else (0, 0, 0, 0)

        ver, reason = verdict(bmrb_res, ca, cb, pdb_res, n_chains, h, s)
        delta = abs(bmrb_res - pdb_res)
        flag  = {"PASS": "✓", "WARN": "~", "FAIL": "✗"}[ver]
        rsn   = f"  {reason}" if reason else ""
        print(f"{bmrb_id:>6} {pdb_id:>6} | {bmrb_res:>5} {ca:>5} {cb:>5} | {pdb_res:>5} {n_chains:>3} {h:>4} {s:>4} | {delta:>4} | {flag} {ver:<4}{rsn}  [{notes}]")

    print("\n\nAdd PASS and WARN (CB>0) entries to SOLID_STATE_ENTRIES in pipeline.py.")


if __name__ == "__main__":
    main()
