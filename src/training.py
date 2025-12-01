import os, sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
import math
import model.pdbutils as PDB
import model.model as model

def train(pdb_list):
    """
    Train an objective function using interatomic distance distributions
    from a dataset of known 3D structures.
    """
    # Initialize counts
    sum_pair_counts = {bp: [0] * model.max_distance for bp in model.base_pairs}
    sum_reference_counts = [0] * model.max_distance

    # Aggregate counts from all PDB files
    for pdb_file in pdb_list:
        atoms = PDB.extract_c3_atoms(pdb_file)

        reference_counts, pair_counts = model.distance_counts(atoms)
        # Add to global counts
        sum_reference_counts = [a + b for a, b in zip(sum_reference_counts, reference_counts)]
        for bp in model.base_pairs:
            sum_pair_counts[bp] = [a + b for a, b in zip(sum_pair_counts[bp], pair_counts[bp])]

    # Compute scores
    reference_distances_freq = model.frequencies(sum_reference_counts)
    scores = {}
    for bp in model.base_pairs:
        distribution = [0.0] * model.max_distance
        pair_distances_freq = model.frequencies(sum_pair_counts[bp])
        for i in range(model.max_distance):
            pair_freq = pair_distances_freq[i]
            ref_freq = reference_distances_freq[i]
            u = model.score(ref_freq, pair_freq)
            if math.isnan(u) or u > model.maximum_score:
                u = model.maximum_score
            distribution[i] = u
        scores[bp] = distribution
    return scores


def run_train():
    # Collect all PDB files from data/pdbs/native/
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

    # Train the profile
    distributions = train(training_instances)

    # Write the profile
    profile_dir = os.path.join("data", "profiles")
    os.makedirs(profile_dir, exist_ok=True)

    for bp in model.base_pairs:
        output_file = os.path.join(profile_dir, bp + ".txt")
        with open(output_file, "w") as f:
            for i in range(model.max_distance):
                f.write(f"{distributions[bp][i]}\n")
    
    print(f"Profiles saved to {profile_dir}")
