import os, sys
sys.path.append("../")
import numpy as np
import matplotlib.pyplot as plt
from utils.pair import set_pairs

nucleotides = ["A", "U", "G", "C"]
base_pairs = set_pairs(nucleotides)
distance_range = list(range(1, 21))  # precomputed once


def make_single_plot(distribution, distance_range, residue_1, residue_2, output_path):
    plt.figure()
    plt.plot(distance_range, distribution)
    plt.xlabel("distance (Å)")
    plt.ylabel("score")
    plt.title(f"Interaction profile {residue_1}-{residue_2}")
    plt.savefig(output_path)   


def make_plot():
    profile_dir = os.path.join("data", "profiles")
    plot_dir = os.path.join("data", "plots")

    if not os.path.isdir(profile_dir):
        raise FileNotFoundError(f"Profiles folder {profile_dir} not found")

    os.makedirs(plot_dir, exist_ok=True)

    for pair in base_pairs:
        filename = os.path.join(profile_dir, f"{pair}.txt")
        if not os.path.exists(filename):
            continue

        # Fast loading
        distribution = np.loadtxt(filename).tolist()

        output_file = os.path.join(plot_dir, f"{pair}.png")
        make_single_plot(distribution, distance_range, pair[0], pair[1], output_file)
        plt.show()
        plt.close()

    print(f"Plots saved to {plot_dir}")
