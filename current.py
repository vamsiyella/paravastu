
import os
import requests
import tempfile
from collections import defaultdict
import numpy as np
import pandas as pd

# Parsers & structure tools
import pynmrstar  # to parse NMR‐STAR entries :contentReference[oaicite:3]{index=3}
#import nmrstarlib  # optional, for bulk handling / JSON conversion :contentReference[oaicite:4]{index=4}
from Bio.PDB import PDBList, PDBParser
from Bio.PDB.DSSP import DSSP

# Optional: caching
from functools import lru_cache

# Configuration: which BMRB entries to include
# Could be a list of IDs you know are solid‐state, or filtered from metadata
BMRB_IDS = [27661
    # e.g. 17561, 17661, etc.
]

# Optional: which atom types you want stats for
ATOM_TYPES = ['C', 'CA', 'CB', 'N']  # maybe also CO, HN etc.

def fetch_bmrb_star(bmr_id):
    url = f"https://bmrb.io/ftp/pub/bmrb/entry_directories/bmr{bmr_id}/bmr{bmr_id}_3.str"
    resp = requests.get(url, timeout=30)
    resp.raise_for_status()
    return resp.text

def parse_shifts(nmrstar_text):
    entry = pynmrstar.Entry.from_string(nmrstar_text)
    shifts = []

    loops = entry.get_loops_by_category('Atom_chem_shift')

    for loop in loops:
        tags = loop.tags  # column names
        tag_index = {tag: i for i, tag in enumerate(tags)}

        # Required columns
        if not all(t in tag_index for t in ['Comp_ID', 'Atom_ID', 'Val', 'Seq_ID']):
            continue

        for row in loop:
            try:
                comp = row[tag_index['Comp_ID']]
                atom = row[tag_index['Atom_ID']]
                shift_val = float(row[tag_index['Val']])
                seq_id = row[tag_index['Seq_ID']]
            except (ValueError, IndexError, TypeError):
                continue

            shifts.append({
                'residue_name': comp,
                'atom': atom,
                'shift': shift_val,
                'seq_id': seq_id
            })

    return pd.DataFrame(shifts)

def residue_shift_coverage(df_shifts):
    """
    For each residue number, report which atoms have shifts.
    """
    coverage = defaultdict(set)

    for _, row in df_shifts.iterrows():
        seq_id = row['seq_id']
        atom = row['atom']
        if atom in ATOM_TYPES:
            coverage[int(seq_id)].add(atom)

    rows = []
    for resnum in sorted(coverage):
        atoms = sorted(coverage[resnum])
        rows.append({
            'seq_id': resnum,
            'atoms_observed': atoms,
            'n_atoms': len(atoms)
        })

    return pd.DataFrame(rows)



@lru_cache()
def fetch_pdb_file(pdb_id, path=None):
    """Download PDB file, return path."""
    pdbl = PDBList()
    if path is None:
        path = tempfile.gettempdir()
    # get_pdb_file will store under PDB structure
    fname = pdbl.retrieve_pdb_file(pdb_id, pdir=path, file_format='pdb')
    return fname

def compute_secondary_structure(pdb_file, chain_id=None):
    """Run DSSP on pdb_file, get mapping of residue (seq_id/resname) → SS class (H, E, etc.)"""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(os.path.basename(pdb_file), pdb_file)
    model = structure[0]  # first model
    # If chain_id not given, take first chain
    if chain_id is None:
        chain = list(model.get_chains())[0]
    else:
        chain = model[chain_id]
    dssp = DSSP(model, pdb_file, dssp='mkdssp')  # requires DSSP executable installed
    # dssp mapping keys are (chain_id, resseq)
    ss_map = {}
    for key, dssp_data in dssp.property_dict.items():
        # key: (chain_id, resseq, insertion code)
        (ch, resnum, icode) = key
        ss = dssp_data['secondary_structure']  # one of H, B, E, G, I, T, S, etc.
        # Optionally map to coarse categories: helix (H, G, I), strand (E, B), coil/other
        if ss in ('H', 'G', 'I'):
            ss_class = 'helix'
        elif ss in ('E', 'B'):
            ss_class = 'strand'
        else:
            ss_class = 'coil'
        # store
        ss_map[(ch, str(resnum))] = ss_class
    return ss_map

def extract_dssp_residue_map(pdb_file, chain_id=None):
    """
    Run DSSP and return a mapping:
    residue_number (int) -> ss_class ('helix', 'strand', 'coil')
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("struct", pdb_file)
    model = structure[0]

    if chain_id is None:
        chain = list(model.get_chains())[0]
        chain_id = chain.id

    dssp = DSSP(model, pdb_file, dssp='mkdssp')

    ss_map = {}

    for (ch, resseq, icode), dssp_data in dssp.property_dict.items():
        if ch != chain_id:
            continue

        ss = dssp_data['secondary_structure']

        if ss in ('H', 'G', 'I'):
            ss_class = 'helix'
        elif ss in ('E', 'B'):
            ss_class = 'strand'
        else:
            ss_class = 'coil'

        ss_map[int(resseq)] = ss_class

    return ss_map

def build_dssp_segments(ss_map):
    """
    Convert per-residue DSSP assignments into contiguous segments.
    Returns a list of dicts:
    {start, end, ss_class}
    """
    segments = []

    sorted_residues = sorted(ss_map.keys())
    if not sorted_residues:
        return segments

    current_start = sorted_residues[0]
    current_ss = ss_map[current_start]
    prev_res = current_start

    for res in sorted_residues[1:]:
        if res == prev_res + 1 and ss_map[res] == current_ss:
            prev_res = res
            continue

        segments.append({
            'start': current_start,
            'end': prev_res,
            'ss_class': current_ss
        })

        current_start = res
        current_ss = ss_map[res]
        prev_res = res

    segments.append({
        'start': current_start,
        'end': prev_res,
        'ss_class': current_ss
    })

    return segments

def get_pdb_id_from_bmrb_text(nmrstar_text):
    """
    Extract PDB ID from an already-loaded BMRB NMR-STAR text.
    Uses frame_list for maximum pynmrstar compatibility.
    """
    entry = pynmrstar.Entry.from_string(nmrstar_text)

    for sf in entry.frame_list:
        if sf.category in ('entry_information', 'deposited_coordinates'):
            try:
                pdb_ids = sf.get_tag('PDB_ID')
                if pdb_ids:
                    return pdb_ids[0]
            except KeyError:
                continue

    return None

def get_uniprot_from_bmrb_text(nmrstar_text):
    """
    Extract UniProt accession from BMRB entry if available.
    """
    entry = pynmrstar.Entry.from_string(nmrstar_text)

    for sf in entry.frame_list:
        if sf.category == 'entry_information':
            try:
                accessions = sf.get_tag('Database_accession_code')
                if accessions:
                    return accessions[0]
            except KeyError:
                continue

    return None






def residue_three_to_one(res3):
    """Convert three‐letter code to one‐letter code, optionally."""
    mapping = {
        'ALA':'A','CYS':'C','ASP':'D','GLU':'E','PHE':'F','GLY':'G',
        'HIS':'H','ILE':'I','LYS':'K','LEU':'L','MET':'M','ASN':'N',
        'PRO':'P','GLN':'Q','ARG':'R','SER':'S','THR':'T','VAL':'V',
        'TRP':'W','TYR':'Y'
    }
    return mapping.get(res3.upper(), 'X')



def extract_sequence_from_bmrb(nmrstar_text):
    """
    Robust sequence extraction by inferring residues
    directly from Atom_chem_shift loops.
    Works for solid-state BMRB entries.
    """
    entry = pynmrstar.Entry.from_string(nmrstar_text)

    seq_records = {}

    # Pull ALL Atom_chem_shift loops, regardless of saveframe
    loops = entry.get_loops_by_category('Atom_chem_shift')

    for loop in loops:
        tags = loop.tags

        if 'Comp_ID' not in tags or 'Seq_ID' not in tags:
            continue

        comp_idx = tags.index('Comp_ID')
        seq_idx = tags.index('Seq_ID')

        for row in loop:
            try:
                resnum = int(row[seq_idx])
                resname = row[comp_idx]
            except (TypeError, ValueError, IndexError):
                continue

            # Only record first occurrence per residue
            if resnum not in seq_records:
                seq_records[resnum] = resname

    if not seq_records:
        return None

    sequence = ''.join(
        residue_three_to_one(seq_records[i])
        for i in sorted(seq_records)
    )

    return sequence

def build_residue_atom_coverage(df_shifts):
    """
    Summarize which atoms are observed for each residue.
    """
    coverage = (
        df_shifts
        .groupby('seq_id')['atom']
        .apply(lambda x: sorted(set(x)))
        .reset_index()
    )
    coverage['n_atoms'] = coverage['atom'].apply(len)
    coverage.rename(columns={'atom': 'atoms_observed'}, inplace=True)
    return coverage


def aggregate_statistics(bmrb_ids):
    """Main function to build statistics table: shift distributions by (residue, atom, secondary structure)."""
    # Data structure: dict keyed by (residue_one, atom, ss_class) -> list of shifts
    data = defaultdict(list)

    for bmr_id in bmrb_ids:
        print(f"\n=== Processing BMRB {bmr_id} ===CHANGE 1")
        print(f"Processing BMRB {bmr_id}...")
        try:
            text = fetch_bmrb_star(bmr_id)
        except Exception as e:
            print(f"  Error fetching {bmr_id}: {e}")
            continue
        df_shifts = parse_shifts(text)
        print(f"Parsed shifts: {len(df_shifts) } CHANGE 2")

        if df_shifts.empty:
            continue

        pdb_id = None

        ss_map = None
        if pdb_id:
            print(f"PDB ID found: {pdb_id} CHANGE 3")

            try:
                pdb_file = fetch_pdb_file(pdb_id)
                ss_map = compute_secondary_structure(pdb_file)
            except Exception as e:
                print(f"  Could not process PDB {pdb_id} secondary structure: {e}")
                ss_map = None

        for idx, row in df_shifts.iterrows():
            res3 = row['residue_name']
            atom = row['atom']
            shift = row['shift']
            seq_id = row['seq_id']
            res1 = residue_three_to_one(res3)
            if atom not in ATOM_TYPES:
                continue
            """if ss_map is not None:
                # need chain id; assume chain is known or single chain
                # key in ss_map: (chain_id, seq_id) -- seq_id matching row['Seq_ID']
                # This matching may need adjusting depending on numbering/insertion codes etc.
                # For simplicity assume chain 'A' or first chain
                key = ('A', str(seq_id))
                ss_class = ss_map.get(key, 'unknown')
            else:
                ss_class = 'unknown' """
            
            ss_class = 'unknown'  # CHANGE 4: skip SS classification for now

            data[(res1, atom, ss_class)].append(shift)

    # Now compute statistics
    rows = []
    for (res1, atom, ss_class), shifts in data.items():
        if len(shifts) == 0:
            continue
        arr = np.array(shifts)
        rows.append({
            'residue': res1,
            'atom': atom,
            'ss_class': ss_class,
            'count': len(arr),
            'mean': np.mean(arr),
            'median': np.median(arr),
            'std': np.std(arr),
            'min': np.min(arr),
            'max': np.max(arr)
        })

    stats_df = pd.DataFrame(rows)
    return stats_df

if __name__ == "__main__":
    test_id = 17561
    text = fetch_bmrb_star(test_id)

    seq = extract_sequence_from_bmrb(text)

    if seq is None:
        print("No sequence found!")
    else:
        print(f"BMRB {test_id} sequence length:", len(seq))
        print(seq)

    df_shifts = parse_shifts(text)

    # Analyze residue-level coverage
    coverage_df = residue_shift_coverage(df_shifts)

    print("\nResidue-level atom coverage (first 10 residues):")
    print(coverage_df.head(10))

    print("\nAverage atoms per residue:")
    print(coverage_df['n_atoms'].mean())    

    df_shifts = parse_shifts(text)
    coverage = build_residue_atom_coverage(df_shifts)

    print("\nResidue-level atom coverage (first 10 residues):")
    print(coverage.head(10))

    print("\nAverage atoms per residue:")
    print(coverage['n_atoms'].mean())

    print("\nFetching PDB ID from BMRB...")
    pdb_id = get_pdb_id_from_bmrb_text(text)

    print("PDB ID:", pdb_id)

    if pdb_id is None:
        print("No PDB linked to this BMRB entry.")
    else:
        pdb_file = fetch_pdb_file(pdb_id)
        ss_map = extract_dssp_residue_map(pdb_file)

        print("\nDSSP residue-level assignments (first 10):")
        for res in sorted(ss_map)[:10]:
            print(f"Residue {res}: {ss_map[res]}")

        segments = build_dssp_segments(ss_map)

        print("\nDSSP segments:")
        for seg in segments[:10]:
            print(seg)
    if pdb_id is None:
        print("\nTrying to extract UniProt ID...")
        uniprot_id = get_uniprot_from_bmrb_text(text)
        print("UniProt ID:", uniprot_id)
