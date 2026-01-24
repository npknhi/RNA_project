# src/plotting.py

import os
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from utils.pair import set_pairs


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

def plot_kde(distance_file="data/distances.csv", bandwidth=0.5):
    """
    Read a TSV/CSV of base pair distances and plot KDEs per base pair.

    Parameters
    ----------
    distance_file : str
        Path to the TSV/CSV file containing base pair distances.
    bandwidth : float
        KDE bandwidth (sigma).
    """

    plot_dir = os.path.join("data", "plots")
    os.makedirs(plot_dir, exist_ok=True)

    df = pd.read_csv(distance_file, sep="\t")

    if df.empty:
        raise RuntimeError(f"No data found in {distance_file}")

    base_pairs = df["base_pair"].unique()
    max_distance = df["distance"].max()
    grid_step = 0.1
    grid = np.arange(0.0, max_distance + grid_step, grid_step)

    plt.figure(figsize=(10, 6))

    for bp in base_pairs:
        distances = df.loc[df["base_pair"] == bp, "distance"].values
        if len(distances) == 0:
            continue

        # Gaussian KDE
        diff = grid[:, None] - distances[None, :]
        kernel = np.exp(-0.5 * (diff / bandwidth) ** 2)
        density = kernel.sum(axis=1)
        density /= len(distances) * bandwidth * np.sqrt(2 * np.pi)

        plt.plot(grid, density, label=bp)

    plt.xlabel("Distance (Å)")
    plt.ylabel("Density")
    plt.title(f"KDE of inter-residue distances from {distance_file}")
    plt.legend()
    plt.tight_layout()

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
