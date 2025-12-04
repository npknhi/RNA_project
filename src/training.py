import os
import argparse
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
    reference_freq = model.frequencies(sum_reference_counts)

    # Compute pair-based distributions
    scores = {}
    for bp in model.base_pairs:
        pair_freq = model.frequencies(sum_pair_counts[bp])

        ref = reference_freq
        pf = pair_freq
        u = np.empty_like(ref)

        mask_zero = (ref == 0) | (pf == 0)

        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = np.divide(pf, ref, out=np.zeros_like(pf), where=~mask_zero)
            u_raw = -np.log(ratio)

        u[:] = u_raw
        u[mask_zero] = model.maximum_score
        u[np.isnan(u) | (u > model.maximum_score)] = model.maximum_score

        scores[bp] = u.tolist()

    return scores


def run_train(train_dir, profile_dir):
    """
    Unified training function using user-specified train dataset/profile output folders.
    """

    if not os.path.isdir(train_dir):
        raise FileNotFoundError(f"Dataset folder {train_dir} not found")

    train_files = [
        os.path.join(train_dir, f)
        for f in os.listdir(train_dir)
        if f.endswith(".pdb")
    ]

    if not train_files:
        raise RuntimeError(f"No PDB files found in {train_dir}")

    # Train
    distributions = train(train_files)

    # Save profile output
    os.makedirs(profile_dir, exist_ok=True)

    for bp in model.base_pairs:
        output_file = os.path.join(profile_dir, f"{bp}.txt")
        with open(output_file, "w") as f:
            for value in distributions[bp]:
                f.write(f"{value}\n")

    print(f"Profiles saved to {profile_dir}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Training module for RNA scoring model")

    parser.add_argument(
        "--trainset",
        default="data/pdbs/train",
        help="Folder containing training PDB structures (default: data/pdbs/train)"
    )

    parser.add_argument(
        "--output",
        default="data/profiles",
        help="Folder where trained model profiles will be written (default: data/profiles)"
    )

    parser.add_argument("--max-distance", type=int, default=None,
                        help="Set maximum allowed distance cutoff (default: 20)")

    parser.add_argument("--position-skip", type=int, default=None,
                        help="Minimum residue separation (default: 4)")

    parser.add_argument("--maximum-score", type=int, default=None,
                        help="Maximum allowed score value in the statistical potential (default: 10)")

    args = parser.parse_args()

    # Apply overrides only if provided
    if args.max_distance is not None:
        model.max_distance = args.max_distance
        model.max_distance_sq = args.max_distance ** 2  

    if args.position_skip is not None:
        model.position_skip = args.position_skip

    if args.maximum_score is not None:
        model.maximum_score = args.maximum_score

    print("Training parameters in use:")
    print("  trainset_dir   =", args.trainset)
    print("  output_dir     =", args.output)
    print("  max_distance   =", model.max_distance)
    print("  position_skip  =", model.position_skip)
    print("  maximum_score  =", model.maximum_score)


    run_train(args.trainset, args.output)
