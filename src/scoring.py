import math
import os, sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
import csv

import model.pair as pair
import model.pdbutils as PDB
import model.model as model
import model.interpolation as interpolation


def score(atoms, reference_distributions):
    """
    Compute the estimated Gibbs free energy of an RNA conformation
    by summing scoring values for all residue distances.
    """
    distances = model.residue_distances(atoms)
    s = 0.0
    for residue_i, residue_j, distance in distances:
        reference_distribution = reference_distributions[pair.normalize_pair(residue_i, residue_j)]
        x0 = math.ceil(distance)
        x1 = math.floor(distance)
        if x0 >= len(reference_distribution):
            x0 = len(reference_distribution) - 1
        if x1 >= len(reference_distribution):
            x1 = len(reference_distribution) - 1
        y0 = reference_distribution[x0]
        y1 = reference_distribution[x1]
        u = interpolation.linear_interpolation(x0, y0, x1, y1, distance)
        s += u
    return s

def run_score():
    profile_dir = os.path.join("data", "profiles")
    if not os.path.isdir(profile_dir):
        raise FileNotFoundError(f"Profiles folder {profile_dir} not found")

    # Load reference distributions
    reference_distributions = {}
    for bp in model.base_pairs:
        filename = os.path.join(profile_dir, bp + ".txt")
        with open(filename, "r") as f:
            reference_distributions[bp] = [float(line.strip()) for line in f if line.strip()]

    # Collect all decoy PDB files
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

    # Compute scores for all decoys
    results = []
    for pdb_file in decoy_files:
        atoms = PDB.extract_c3_atoms(pdb_file)
        s = score(atoms, reference_distributions)
        print(f" - {os.path.basename(pdb_file)}: {s:.4f}")
        results.append((os.path.basename(pdb_file), s))


    # Save results to CSV
    results_dir = os.path.join("data", "scores")
    os.makedirs(results_dir, exist_ok=True)
    output_file = os.path.join(results_dir, "scores.csv")

    with open(output_file, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["pdb_file", "score"])
        writer.writerows(results)
    print(f"Scores saved to {output_file}")


