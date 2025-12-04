import math
import os
import argparse
import datetime
import csv
import numpy as np

import utils.pair as pair
import utils.pdb as pdb
import utils.model as model
import utils.interpolation as interpolation


def score(atoms, reference_distributions):
    """
    Compute the estimated Gibbs free energy of an RNA conformation.
    """
    distances = model.residue_distances(atoms)
    s = 0.0

    norm_pair = pair.normalize_pair
    lin_interp = interpolation.linear_interpolation

    for residue_i, residue_j, d in distances:
        key = norm_pair(residue_i, residue_j)
        rd = reference_distributions[key]

        x1 = int(math.floor(d))
        x0 = x1 + 1

        if x1 >= len(rd):
            x1 = len(rd) - 1
        if x0 >= len(rd):
            x0 = len(rd) - 1

        y1 = rd[x1]
        y0 = rd[x0]

        s += lin_interp(x0, y0, x1, y1, d)

    return s


def run_score(model_dir, testset_dir, output_dir):

    # === Load reference profiles ===
    if not os.path.isdir(model_dir):
        raise FileNotFoundError(f"Model folder {model_dir} not found")

    reference_distributions = {}
    for bp in model.base_pairs:
        filename = os.path.join(model_dir, f"{bp}.txt")
        reference_distributions[bp] = np.loadtxt(filename).tolist()


    # === Load test PDBs ===
    if not os.path.isdir(testset_dir):
        raise FileNotFoundError(f"Test dataset folder {testset_dir} not found")

    test_files = [
        os.path.join(testset_dir, f)
        for f in os.listdir(testset_dir)
        if f.endswith(".pdb")
    ]

    if not test_files:
        raise RuntimeError(f"No PDB files found in {testset_dir}")


    # === Score all test structures ===
    results = []
    for pdb_file in test_files:
        atoms = pdb.extract_c3_atoms(pdb_file)
        s = score(atoms, reference_distributions)
        print(f" - {os.path.basename(pdb_file)}: {s:.4f}")
        results.append((os.path.basename(pdb_file), s))


    # === Save results ===
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(
        output_dir,
        f"{datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')}_scores.csv"
    )

    with open(output_file, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["pdb_file", "score"])
        writer.writerows(results)

    print(f"Scores saved to {output_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Scoring module for RNA structures")

    parser.add_argument(
        "--model",
        default="data/profiles",
        help="Folder containing trained model profiles (default: data/profiles)"
    )

    parser.add_argument(
        "--testset",
        default="data/pdbs/test",
        help="Folder containing test PDB structures (default: data/pdbs/test)"
    )

    parser.add_argument(
        "--output",
        default="data/scores",
        help="Folder to store scoring results (default: data/scores)"
    )

    args = parser.parse_args()

    print("Scoring parameters:")
    print("  model_dir     =", args.model)
    print("  testset_dir   =", args.testset)
    print("  output_dir    =", args.output)

    run_score(args.model, args.testset, args.output)
