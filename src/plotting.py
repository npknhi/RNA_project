import os
import numpy as np
import matplotlib.pyplot as plt
from utils.pair import set_pairs

nucleotides = ["A", "U", "G", "C"]
base_pairs = set_pairs(nucleotides)
distance_range = list(range(1, 21))


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

        # Load distribution
        distribution = np.loadtxt(filename).tolist()

        plt.figure(figsize=(6, 4))
        plt.plot(distance_range, distribution)
        plt.title(f"Interaction profile {pair}", fontsize=12)
        plt.xlabel("Distance (Å)")
        plt.ylabel("Score")
        plt.grid(True)

        output_file = os.path.join(plot_dir, f"{pair}.png")
        plt.savefig(output_file, dpi=300)
        plt.close()

    print(f"\nAll plots saved to {plot_dir}")


if __name__ == "__main__":
    make_plot()
