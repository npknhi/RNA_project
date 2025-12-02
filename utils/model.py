import math
import numpy as np
from utils.pair import set_pairs, normalize_pair

# Parameters
nucleotides = ["A", "U", "G", "C"]
base_pairs = set_pairs(nucleotides)
max_distance = 20
max_distance_sq = max_distance * max_distance
position_skip = 4
maximum_score = 10


def residue_distances(atoms):
    """
    Compute residue-residue distances using NumPy vectorization.
    Returns list of (res_i, res_j, distance)
    """
    atoms = np.array(atoms, dtype=object)

    chains = atoms[:, 0]
    residues = atoms[:, 1]
    coords = atoms[:, 2:].astype(float)

    n = len(atoms)
    distances = []

    for i in range(n - position_skip):
        chain_i = chains[i]
        res_i = residues[i]
        coord_i = coords[i]

        j_range = np.arange(i + position_skip, n)
        same_chain = (chains[j_range] == chain_i)
        j_valid = j_range[same_chain]

        if j_valid.size == 0:
            continue

        coords_j = coords[j_valid]
        diff = coords_j - coord_i
        d2 = np.sum(diff * diff, axis=1)

        mask = d2 < max_distance_sq
        if not np.any(mask):
            continue

        j_kept = j_valid[mask]
        d2_kept = d2[mask]

        for j_idx, d2_val in zip(j_kept, d2_kept):
            distances.append((res_i, residues[j_idx], math.sqrt(d2_val)))

    return distances


def distance_counts(atoms):
    """
    Compute counts of atomic distances.
    """
    pair_counts = {bp: [0] * max_distance for bp in base_pairs}
    reference_counts = [0] * max_distance

    for res_i, res_j, distance in residue_distances(atoms):
        bin_index = int(distance)
        base_pair = normalize_pair(res_i, res_j)

        reference_counts[bin_index] += 1
        pair_counts[base_pair][bin_index] += 1

    return reference_counts, pair_counts


def frequencies(counts):
    """
    Normalize counts to frequencies.
    """
    counts = np.asarray(counts, dtype=float)
    total = counts.sum()
    if total == 0:
        return np.zeros_like(counts)
    return counts / total

def score(reference_frequency, pair_frequency):
    """
    Compute u_ij = -log(f_ij / f_xx)
    """
    if reference_frequency == 0 or pair_frequency == 0:
        return float("inf")
    return -math.log(pair_frequency / reference_frequency)
