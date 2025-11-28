import sys, os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import numpy as np
import matplotlib
matplotlib.use("Agg")   # Fully headless (no display required)
import matplotlib.pyplot as plt

PAIRS = ["AA","AU","AC","AG","UU","UC","UG","CC","CG","GG"]


def plot_profiles(rna_name, base_folder="data/profiles", out_root="plot"):
    """
    Plot scoring profiles for a single RNA target.
    """

    profile_folder = os.path.join(base_folder, rna_name)
    if not os.path.isdir(profile_folder):
        raise FileNotFoundError(f"Profile folder does not exist: {profile_folder}")

    out_folder = os.path.join(out_root, rna_name)
    os.makedirs(out_folder, exist_ok=True)

    print(f"Plotting profiles for RNA: {rna_name}")
    print(f"Input:  {profile_folder}")
    print(f"Output: {out_folder}")

    for p in PAIRS:
        filepath = os.path.join(profile_folder, f"{p}.txt")

        if not os.path.isfile(filepath):
            print(f"[WARN] Missing profile file: {filepath}")
            continue

        profile = np.loadtxt(filepath)

        plt.figure(figsize=(6, 4))
        plt.plot(range(len(profile)), profile, label=p, linewidth=2)
        plt.title(f"Scoring Profile: {p}")
        plt.xlabel("Distance bin (Å)")
        plt.ylabel("Score")
        plt.grid(True)

        save_path = os.path.join(out_folder, f"{p}.png")
        plt.savefig(save_path, dpi=200, bbox_inches='tight')
        plt.close()

        print(f"✔ Saved: {save_path}")

    fig, axes = plt.subplots(5, 2, figsize=(10, 14))
    axes = axes.flatten()

    for i, p in enumerate(PAIRS):
        filepath = os.path.join(out_folder, f"{p}.png")

        if not os.path.isfile(filepath):
            axes[i].set_title(f"{p} (missing)")
            axes[i].axis("off")
            continue

        img = plt.imread(filepath)
        axes[i].imshow(img)
        axes[i].set_title(p)
        axes[i].axis("off")

    plt.tight_layout()

    combined_path = os.path.join(out_folder, f"{rna_name}_combined.png")
    plt.savefig(combined_path, dpi=200, bbox_inches='tight')
    plt.close()

    print(f"✔ Combined figure saved: {combined_path}")


def main():
    if len(sys.argv) < 2:
        print("Usage: python plotting.py <rna_name>")
        sys.exit(1)

    rna_name = sys.argv[1]
    plot_profiles(rna_name)


if __name__ == "__main__":
    main()
