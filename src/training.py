import os
import argparse
import numpy as np
from math import ceil

import utils.rna_extractor as rna_extractor
import utils.model as model

def train(struct_list):
    """
    Train an objective function by computing interatomic distance distributions from a set of RNA structures.

    Parameters
    ----------
    struct_list : list
        List of file paths to PDB/CIF structures used for training.

    Returns
    -------
    dict
        Dictionary of score distributions for each base pair.
    """

    num_bins = model.num_bins

    # Initialize global counts 
    sum_reference_counts = np.zeros(num_bins, dtype=float)     # >>> changed
    sum_pair_counts = {bp: np.zeros(num_bins, dtype=float) for bp in model.base_pairs}  # >>> changed

    # Aggregate counts from all PDB/CIF files
    for struct_file in struct_list:
        atoms = rna_extractor.extract_c3_atoms(struct_file)

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
    Train score distributions from a dataset of PDB structures and save the results to profile files.

    Parameters
    ----------
    train_dir : str
        Path to the folder containing PDB files used for training.
    profile_dir : str
        Path to the folder where the computed profile `.txt` files will be saved.

    Returns
    -------
    None
    """

    if not os.path.isdir(train_dir):
        raise FileNotFoundError(f"Dataset folder {train_dir} not found")

    train_files = [
        os.path.join(train_dir, f)
        for f in os.listdir(train_dir)
        if f.lower().endswith((".pdb", ".cif", ".mmcif"))
    ]

    if not train_files:
        raise RuntimeError(f"No PDB/CIF files found in {train_dir}")

    # Train
    distributions = train(train_files)

    # Save profile output
    os.makedirs(profile_dir, exist_ok=True)

    # for bp in model.base_pairs:
    #     output_file = os.path.join(profile_dir, f"{bp}.txt")
    #     with open(output_file, "w") as f:
    #         for value in distributions[bp]:
    #             f.write(f"{value}\n")

    for bp in model.base_pairs:
        output_file = os.path.join(profile_dir, f"{bp}.txt")
        with open(output_file, "w") as f:

            # >>> added: compute bin centers
            distance_range = (np.arange(model.num_bins) + 0.5) * model.bin_width

            for dist, value in zip(distance_range, distributions[bp]):
                f.write(f"{dist:.6f}\t{value}\n")

    print(f"Profiles saved to {profile_dir}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Training module for RNA scoring model")

    parser.add_argument(
        "--trainset",
        default="data/structures/train",
        help="Folder containing training PDB/CIF structures (default: data/structures/train)"
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

    parser.add_argument("--bin-width", type=float, default=1.0,
                        help="Histogram bin width for distance distributions (default: 1.0 Å)")


    args = parser.parse_args()

    # Apply overrides only if provided
    if args.max_distance is not None:
        model.max_distance = args.max_distance
        model.max_distance_sq = args.max_distance ** 2  

    if args.position_skip is not None:
        model.position_skip = args.position_skip

    if args.maximum_score is not None:
        model.maximum_score = args.maximum_score

    if args.bin_width is not None:
        model.bin_width = args.bin_width
        model.num_bins = ceil(model.max_distance / model.bin_width)

    print("Training parameters in use:")
    print("  trainset_dir   =", args.trainset)
    print("  output_dir     =", args.output)
    print("  max_distance   =", model.max_distance)
    print("  position_skip  =", model.position_skip)
    print("  maximum_score  =", model.maximum_score)
    print("  bin_width      =", model.bin_width)                     
    print("  num_bins       =", model.num_bins)                      

    run_train(args.trainset, args.output)
