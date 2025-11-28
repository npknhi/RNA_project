import sys, os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from model.pdbutils import extract_c3_atoms
from model.pair import normalize_pair
import numpy as np

# Constants
MAX_DIST = 20          # maximum distance considered (Å)
BIN_SIZE = 1           # 1 Å per bin
N_BINS = MAX_DIST      # 20 bins from 0 to 20 Å
EPS = 1e-12            # small epsilon to avoid division by zero

# 10 RNA base-pair categories used for scoring
PAIRS = [
    "AA","AU","AC","AG",
    "UU","UC","UG",
    "CC","CG","GG"
]

def init_counts():
    """
    Initialize distance histograms and counters for:
        - each base pair (AA, AU, AC, ...)
        - reference distribution (all pairs pooled)
    """
    return (
        {p: np.zeros(N_BINS) for p in PAIRS},   # pair histograms
        {p: 0 for p in PAIRS},                  # total counts per pair
        np.zeros(N_BINS),                       # reference histogram
        0                                       # total reference counts
    )

def process_pdb_file(path, pair_counts, pair_total, ref_counts, ref_total):
    """
    Extract C3' atoms, compute all valid pairwise distances (i,i+4),
    update pair histograms and reference histogram.
    """
    atoms = extract_c3_atoms(path)
    n = len(atoms)

    for i in range(n):
        for j in range(i + 4, n):  # residues must be at least 4 positions apart
            b1, c1 = atoms[i]
            b2, c2 = atoms[j]

            dist = np.linalg.norm(c1 - c2)
            if dist >= MAX_DIST:
                continue

            bin_index = int(dist // BIN_SIZE)
            pair = normalize_pair(b1, b2)

            # Update pair histogram
            if pair in pair_counts:
                pair_counts[pair][bin_index] += 1
                pair_total[pair] += 1

            # Update reference histogram
            ref_counts[bin_index] += 1
            ref_total += 1

    return pair_counts, pair_total, ref_counts, ref_total

def compute_profiles(pair_counts, pair_total, ref_counts, ref_total, outfolder):
    """
    Compute scoring profiles as:
        score = -log( observed(pair) / reference )
    Handle NaN/inf and clip to range [0,10].
    """

    ref_freq = ref_counts / (ref_total + EPS)

    for p in PAIRS:
        # Observed frequencies
        obs_freq = np.divide(
            pair_counts[p],
            pair_total[p] if pair_total[p] > 0 else 1,
            out=np.zeros_like(pair_counts[p]),
            where=True
        )

        # Safe ratio
        ratio = np.divide(
            obs_freq,
            ref_freq,
            out=np.zeros_like(obs_freq),
            where=(ref_freq != 0)
        )

        profile = -np.log(ratio + EPS)
        profile = np.nan_to_num(profile, nan=10.0, posinf=10.0, neginf=10.0)
        profile = np.clip(profile, 0, 10)

        # Save profile to file
        np.savetxt(os.path.join(outfolder, f"{p}.txt"), profile)

def main(rna_name):
    """
    Train statistical potential for a single RNA target.

    Input:
        rna_name = e.g. "5k7c"

    Reads:
        data/pdbs/<rna_name>/native.pdb

    Writes:
        data/profiles/<rna_name>/AA.txt ... GG.txt
    """

    base_folder = "data/pdbs"
    rna_folder = os.path.join(base_folder, rna_name)

    native_path = os.path.join(rna_folder, "native.pdb")
    if not os.path.isfile(native_path):
        raise FileNotFoundError(f"Native PDB not found: {native_path}")

    outfolder = os.path.join("data/profiles", rna_name)
    os.makedirs(outfolder, exist_ok=True)

    pair_counts, pair_total, ref_counts, ref_total = init_counts()

    print(f"Processing native: {native_path}")
    pair_counts, pair_total, ref_counts, ref_total = process_pdb_file(
        native_path, pair_counts, pair_total, ref_counts, ref_total
    )

    compute_profiles(pair_counts, pair_total, ref_counts, ref_total, outfolder)
    print(f"✔ Training completed for {rna_name}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python training.py <rna_name>")
        sys.exit(1)

    rna_name = sys.argv[1]
    main(rna_name)
