import sys, os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from model.pdbutils import extract_c3_atoms
from model.model import RNAModel
import numpy as np
import csv

MAX_DIST = 20

def score_structure(pdb_path, model):
    """
    Compute the pseudo-energy score of a single RNA PDB file.
    """
    atoms = extract_c3_atoms(pdb_path)
    n = len(atoms)
    total = 0.0

    for i in range(n):
        for j in range(i + 4, n):  # must be >= 4 positions apart
            b1, c1 = atoms[i]
            b2, c2 = atoms[j]

            d = np.linalg.norm(c1 - c2)
            if d >= MAX_DIST:
                continue

            total += model.score_distance(b1, b2, d)

    return total

def main(rna_name):
    """
    Evaluate ALL non-native PDB models of a given RNA.
    Read:
        data/pdbs/<rna_name>/
            - native.pdb
            - model_*.pdb
    Write:
        data/results/<rna_name>.csv
    Print:
        Table of model → score
    """
    # folders
    base = "data/pdbs"
    rna_folder = os.path.join(base, rna_name)

    if not os.path.isdir(rna_folder):
        raise FileNotFoundError(f"RNA folder not found: {rna_folder}")

    # load trained profiles
    model = RNAModel(os.path.join("data/profiles", rna_name))

    # list models (exclude native.pdb)
    pdb_files = sorted(
        f for f in os.listdir(rna_folder)
        if f.endswith(".pdb") and f != "native.pdb"
    )

    if not pdb_files:
        print("No model files found.")
        return

    # prepare output result file
    os.makedirs("data/results", exist_ok=True)
    out_csv = os.path.join("data/results", f"{rna_name}.csv")

    # compute scores
    results = []

    for fname in pdb_files:
        path = os.path.join(rna_folder, fname)
        score = score_structure(path, model)
        results.append((fname, score))

    # write csv
    with open(out_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["model", "score"])
        writer.writerows(results)

    # print table to screen
    print(f"\n=== Evaluation Results for {rna_name} ===")
    print(f"{'MODEL':<20} SCORE")
    print("-" * 32)

    for name, sc in results:
        print(f"{name:<20} {sc:.4f}")

    print(f"\nCSV saved to: {out_csv}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python evaluation.py <rna_name>")
        sys.exit(1)

    main(sys.argv[1])
