import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
from utils.pair import set_pairs
import utils.model as model
import utils.rna_extractor as rna_extractor


nucleotides = ["A", "U", "G", "C"]
base_pairs = set_pairs(nucleotides)

def make_plot():
    """
    Generate and save plots of interaction profiles for each base pair.

    Parameters
    ---------------
    None

    Returns
    -------
    None
    """
    profile_dir = os.path.join("data", "profiles")
    plot_dir = os.path.join("data", "plots")

    if not os.path.isdir(profile_dir):
        raise FileNotFoundError(f"Profiles folder {profile_dir} not found")

    os.makedirs(plot_dir, exist_ok=True)

    for pair in base_pairs:
        filename = os.path.join(profile_dir, f"{pair}.txt")
        if not os.path.exists(filename):
            continue

        # Load distribution
        data = np.loadtxt(filename)
        distance = data[:, 0]
        score = data[:, 1]

        plt.figure(figsize=(6, 4))
        plt.plot(distance, score)
        plt.xticks(distance[::2], rotation=45)
        plt.title(f"Interaction profile {pair}", fontsize=12)
        plt.xlabel("Distance (Å)")
        plt.ylabel("Score")
        plt.grid(True)

        output_file = os.path.join(plot_dir, f"{pair}.png")
        plt.savefig(output_file, dpi=300)
        plt.close()

    print(f"\nAll plots saved to {plot_dir}")


def plot_kde(structure_dir="data/structures/test", bandwidth=0.5):
    """
    Compute and display KDEs of inter-residue distances across all RNA structures
    in the given directory.

    Parameters
    ----------
    structure_dir : str
        Folder containing PDB/CIF files to include.
    bandwidth : float
        Bandwidth parameter for KDE smoothing.
    """

    plot_dir = os.path.join("data", "plots")
    os.makedirs(plot_dir, exist_ok=True)
    
    # --- Collect all structure files ---
    if not os.path.isdir(structure_dir):
        raise FileNotFoundError(f"Directory not found: {structure_dir}")

    structure_files = [
        os.path.join(structure_dir, f)
        for f in os.listdir(structure_dir)
        if f.lower().endswith((".pdb", ".cif", ".mmcif"))
    ]

    if not structure_files:
        raise RuntimeError(f"No PDB/CIF files found in {structure_dir}")

    # --- Aggregate distances per base pair ---
    # Initialize a dictionary to store distances for each base pair
    all_distances = {bp: [] for bp in model.base_pairs}

    for struct_file in structure_files:
        atoms = rna_extractor.extract_c3_atoms(struct_file)
        distances = model.residue_distances(atoms)

        for res_i, res_j, d in distances:
            bp = model.normalize_pair(res_i, res_j)
            all_distances[bp].append(d)

    # --- Compute KDEs ---
    kde_data = {}
    grid_step = 0.1
    max_distance = model.max_distance
    grid = np.arange(0.0, max_distance, grid_step)

    for bp, dist_list in all_distances.items():
        if not dist_list:
            kde_data[bp] = np.zeros((len(grid), 2))
            continue

        distances = np.asarray(dist_list)
        diff = grid[:, None] - distances[None, :]
        kernel = np.exp(-0.5 * (diff / bandwidth) ** 2)
        density = kernel.sum(axis=1)
        density /= (len(distances) * bandwidth * (2 * np.pi) ** 0.5)
        kde_data[bp] = np.column_stack((grid, density))

    # --- Plot ---
    plt.figure(figsize=(10, 6))

    colors = {
        "AA": "red", "AU": "orange", "AG": "green", "AC": "blue",
        "CC": "red", "CG": "orange", "CU": "green",
        "GG": "red", "GU": "orange",
        "UU": "red"
    }

    for bp, data in kde_data.items():
        x = data[:, 0]
        y = data[:, 1]
        if y.sum() == 0:
            continue

        plt.plot(x, y, label=bp, color=colors.get(bp, "black"))

    plt.xlabel("Distance (Å)")
    plt.ylabel("Density")
    plt.title(f"KDE of inter-residue distances ({len(structure_files)} structures)")
    plt.legend()
    plt.tight_layout()
    plt.show()

    output_file = os.path.join(plot_dir, "KDE_all.png")
    plt.savefig(output_file, dpi=300)
    plt.close()
    print(f"\n{output_file} saved to {plot_dir}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Scoring module for RNA structures")
    
    parser.add_argument("--make-kde", action="store_true",
                        help="Display KDE plots (default: False)")
    args = parser.parse_args()
    run_kde = args.make_kde
    
    print("KDE plotting:  ", run_kde)
   
    make_plot()

    if run_kde:
        plot_kde()
