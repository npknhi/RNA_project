import math
import os, sys
sys.path.append("../")
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

    # Local binding for speed
    norm_pair = pair.normalize_pair
    lin_interp = interpolation.linear_interpolation
    ref = reference_distributions

    for residue_i, residue_j, d in distances:
        key = norm_pair(residue_i, residue_j)
        rd = ref[key]

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


def run_score():
    profile_dir = os.path.join("data", "profiles")
    if not os.path.isdir(profile_dir):
        raise FileNotFoundError(f"Profiles folder {profile_dir} not found")

    reference_distributions = {}
    for bp in model.base_pairs:
        filename = os.path.join(profile_dir, bp + ".txt")
        reference_distributions[bp] = np.loadtxt(filename).tolist()

    decoy_dir = os.path.join("data", "pdbs", "decoys")
    if not os.path.isdir(decoy_dir):
        raise FileNotFoundError(f"Decoy folder {decoy_dir} not found")

    decoy_files = [
        os.path.join(decoy_dir, f)
        for f in os.listdir(decoy_dir)
        if f.endswith(".pdb")
    ]

    if not decoy_files:
        raise RuntimeError(f"No PDB files found in {decoy_dir}")

    results = []
    for pdb_file in decoy_files:
        atoms = pdb.extract_c3_atoms(pdb_file)
        s = score(atoms, reference_distributions)
        print(f" - {os.path.basename(pdb_file)}: {s:.4f}")
        results.append((os.path.basename(pdb_file), s))

    # Save results
    results_dir = os.path.join("data", "scores")
    os.makedirs(results_dir, exist_ok=True)
    output_file = os.path.join(results_dir, "scores.csv")

    with open(output_file, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["pdb_file", "score"])
        writer.writerows(results)

    print(f"Scores saved to {output_file}")
