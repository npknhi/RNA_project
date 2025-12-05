import math
import os
import argparse
import datetime
import csv
import numpy as np

import utils.pair as pair
import utils.rna_extractor as rna_extractor
import utils.model as model
import utils.interpolation as interpolation

def score(atoms, reference_distributions):
    """
    Compute an estimated Gibbs free energy score for an RNA conformation.

    Parameters
    ----------
    atoms : list
        List of atom objects representing the RNA conformation.
    reference_distributions : dict
        Dictionary mapping normalized residue pairs to reference distance distributions.

    Returns
    -------
    float
        Estimated Gibbs free energy of the RNA conformation.
    """
    distances = model.residue_distances(atoms)
    s = 0.0

    norm_pair = pair.normalize_pair
    lin_interp = interpolation.linear_interpolation

    for residue_i, residue_j, d in distances:
        key = norm_pair(residue_i, residue_j)
        rd = reference_distributions[key]

        centers = rd[:, 0]
        scores  = rd[:, 1]

        # Find where d fits in bin centers
        idx = np.searchsorted(centers, d)

        # Clamp to edges
        if idx == 0:
            s += scores[0]
            continue
        if idx >= len(centers):
            s += scores[-1]
            continue

        x0 = centers[idx - 1]
        y0 = scores[idx - 1]
        x1 = centers[idx]
        y1 = scores[idx]

        s += lin_interp(x0, y0, x1, y1, d)

    return s

def run_score(model_dir, testset_dir, output_dir):
    """
    Score a set of RNA structures against reference profiles and save the results.

    Parameters
    ----------
    model_dir : str
        Path to the folder containing reference profile `.txt` files for base pairs.
    testset_dir : str
        Path to the folder containing PDB/CIF files of test RNA structures.
    output_dir : str
        Path to the folder where scoring results CSV will be saved.
    
    Returns
    -------
    None
    """

    # === Load reference profiles ===
    if not os.path.isdir(model_dir):
        raise FileNotFoundError(f"Model folder {model_dir} not found")

    reference_distributions = {}
    for bp in model.base_pairs:
        filename = os.path.join(model_dir, f"{bp}.txt")
        data = np.loadtxt(filename)
        reference_distributions[bp] = data
        # reference_distributions[bp] = np.loadtxt(filename).tolist()

    # === Load test PDBs/CIFs ===
    if not os.path.isdir(testset_dir):
        raise FileNotFoundError(f"Test dataset folder {testset_dir} not found")

    test_files = [
        os.path.join(testset_dir, f)
        for f in os.listdir(testset_dir)
        if f.lower().endswith((".pdb", ".cif", ".mmcif"))
    ]

    if not test_files:
        raise RuntimeError(f"No PDB/CIF files found in {testset_dir}")

    # === Score all test structures ===
    results = []
    for struct_file in test_files:
        atoms = rna_extractor.extract_c3_atoms(struct_file)
        s = score(atoms, reference_distributions)
        print(f" - {os.path.basename(struct_file)}: {s:.4f}")
        results.append((os.path.basename(struct_file), s))

    # === Save results ===
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(
        output_dir,
        f"{datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')}_scores.csv"
    )

    with open(output_file, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["struct_file", "score"])
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
        default="data/structures/test",
        help="Folder containing test PDB/CIF structures (default: data/structures/test)"
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
