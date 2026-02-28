import pandas as pd
from pathlib import Path

from dssp_module import download_pdb, extract_dssp  

def load_sequence():
    """
    Returns protein sequence for BMRB 17561.
    Later this will be dynamically fetched.
    """
    return "VLDLDVRTCPCGPGGKGRCFGICCGDELGCFVGTAEALRCQEENYLPCQSGQKPCGSGGRCAAAGICCPDGCHEDPACDPEAAFS"


def load_coverage_data():
    """
    Loads coverage CSV from data directory automatically.
    """
    data_path = Path(__file__).resolve().parent.parent / "data"

    # Auto-detect first CSV in data folder
    csv_files = list(data_path.glob("*.csv"))
    if not csv_files:
        raise FileNotFoundError("No CSV file found in data/ directory.")

    df = pd.read_csv(csv_files[0])
    return df


def compute_average_atoms(df):
    """
    Computes average number of unique atom types observed per residue.
    """
    residue_atom_counts = df.groupby("residue")["atom"].nunique()
    return residue_atom_counts.mean()


def summarize_coverage(df):
    print("\nResidue-level atom coverage (first 10 residues):")
    print(df.head(10))

    avg_atoms = compute_average_atoms(df)
    print("\nAverage atoms per residue:")
    print(avg_atoms)


def main():
    sequence = load_sequence()
    print(f"BMRB 17561 sequence length: {len(sequence)}")
    print(sequence)

    df = load_coverage_data()
    summarize_coverage(df)

    
    pdb_id = "2LBH"

    print("\nDownloading PDB...")
    pdb_path = download_pdb(pdb_id)

    print("Extracting DSSP...")
    ss_map = extract_dssp(pdb_path)

    print("\nFirst 10 DSSP assignments:")
    for k in list(ss_map.keys())[:10]:
        print(k, ss_map[k])

if __name__ == "__main__":
    main()
     