import os
import requests
from Bio.PDB import PDBParser, DSSP
from Bio import pairwise2
from Bio.Seq import Seq
from pathlib import Path




def download_pdb(pdb_id):
    project_root = Path(__file__).resolve().parent.parent
    data_dir = project_root / "data"
    data_dir.mkdir(exist_ok=True)

    pdb_path = data_dir / f"{pdb_id}.pdb"

    if pdb_path.exists():
        print("PDB already exists locally.")
        return str(pdb_path)

    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    response = requests.get(url)

    with open(pdb_path, "w") as f:
        f.write(response.text)

    return str(pdb_path)


def extract_dssp(pdb_path):
    """
    Runs DSSP and returns residue -> secondary structure mapping.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_path)
    model = structure[0]

    dssp = DSSP(model, pdb_path)

    ss_map = {}
    for key in dssp.keys():
        chain_id, res_id = key
        residue_number = res_id[1]
        ss = dssp[key][2]
        ss_map[residue_number] = ss

    return ss_map


def align_sequences(seq1, seq2):
    """
    Aligns two sequences and returns alignment.
    """
    alignments = pairwise2.align.globalxx(seq1, seq2)
    return alignments[0]