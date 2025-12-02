import os, sys
sys.path.append("../")
import numpy as np
import utils.pdb as pdb
import utils.model as model


def train(pdb_list):
    """
    Train an objective function using interatomic distance distributions
    from a dataset of known 3D structures.
    """
    max_d = model.max_distance

    # Initialize global counts 
    sum_reference_counts = np.zeros(max_d, dtype=float)
    sum_pair_counts = {bp: np.zeros(max_d, dtype=float) for bp in model.base_pairs}

    # Aggregate counts from all PDB files
    for pdb_file in pdb_list:
        atoms = pdb.extract_c3_atoms(pdb_file)
        reference_counts, pair_counts = model.distance_counts(atoms)

        sum_reference_counts += np.asarray(reference_counts, dtype=float)
        for bp in model.base_pairs:
            sum_pair_counts[bp] += np.asarray(pair_counts[bp], dtype=float)

    # Compute reference frequency distribution
    reference_distances_freq = model.frequencies(sum_reference_counts)

    scores = {}
    for bp in model.base_pairs:
        # Compute the pair-specific frequency distribution
        pair_distances_freq = model.frequencies(sum_pair_counts[bp])

        ref = reference_distances_freq
        pair_f = pair_distances_freq
        u = np.empty_like(ref)

        # Positions where reference or pair frequency is zero
        mask_zero = (ref == 0) | (pair_f == 0)

        # Compute raw scores where valid
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = np.divide(pair_f, ref, out=np.zeros_like(pair_f), where=~mask_zero)
            u_raw = -np.log(ratio)

        u[:] = u_raw

        # Replace invalid positions with maximum score
        u[mask_zero] = model.maximum_score

        # Replace NaN or excessive values with maximum score
        mask_nan_or_big = np.isnan(u) | (u > model.maximum_score)
        u[mask_nan_or_big] = model.maximum_score

        scores[bp] = u.tolist()

    return scores


def run_train():
    # Collect all training PDB files from data/pdbs/native/
    dataset_dir = os.path.join("data", "pdbs", "native")
    if not os.path.isdir(dataset_dir):
        raise FileNotFoundError(f"Dataset folder {dataset_dir} not found")

    training_instances = [
        os.path.join(dataset_dir, f)
        for f in os.listdir(dataset_dir)
        if f.endswith(".pdb")
    ]

    if not training_instances:
        raise RuntimeError(f"No PDB files found in {dataset_dir}")

    # Train statistical profiles
    distributions = train(training_instances)

    # Save trained profiles
    profile_dir = os.path.join("data", "profiles")
    os.makedirs(profile_dir, exist_ok=True)

    for bp in model.base_pairs:
        output_file = os.path.join(profile_dir, bp + ".txt")
        with open(output_file, "w") as f:
            for value in distributions[bp]:
                f.write(f"{value}\n")

    print(f"Profiles saved to {profile_dir}")
