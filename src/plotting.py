import os, sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
import argparse
import matplotlib.pyplot as plt
from model.pair import set_pairs

# Define nucleotides
nucleotides = ["A", "U", "G", "C"]

def make_single_plot(distribution, distance_range, residue_1, residue_2, output_path):
    """
    Plot the scoring profile (score vs interatomic distance).
    """
    plt.figure()
    plt.plot(distance_range, distribution)
    plt.xlabel("distance (Å)")
    plt.ylabel("score")
    plt.title(f"Interaction profile {residue_1}-{residue_2}")
    plt.savefig(output_path)
    plt.show()
    plt.close()


def make_plot():
    """
    Generate plots for the given RNA.
    """
    profile_dir = os.path.join("data", "profiles")
    plot_dir = os.path.join("data", "plots")

    if not os.path.isdir(profile_dir):
        raise FileNotFoundError(f"Profiles folder {profile_dir} not found")

    os.makedirs(plot_dir, exist_ok=True)

    # Loop over all base pairs
    for pair in set_pairs(nucleotides):
        filename = os.path.join(profile_dir, f"{pair}.txt")
        if not os.path.exists(filename):
            continue

        # Read distribution values
        with open(filename, "r") as f:
            distribution = [float(line.strip()) for line in f if line.strip()]

        # Plot and save into the RNA-specific folder
        output_file = os.path.join(plot_dir, f"{pair}.png")
        make_single_plot(distribution, list(range(1, 21)), pair[0], pair[1], output_file)
    
    print(f"Plots saved to {plot_dir}")
