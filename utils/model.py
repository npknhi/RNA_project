# utils/model.py

import math
from math import ceil
import numpy as np
import os
import csv
from utils.pair import set_pairs, normalize_pair
from utils.rna_extractor import extract_c3_atoms, extract_all_atom_residues

# Parameters
nucleotides = ("A", "U", "G", "C")
base_pairs = set_pairs(nucleotides)
max_distance = 20
max_distance_sq = max_distance * max_distance
position_skip = 4
maximum_score = 10
bin_width = 1.0                      
num_bins = ceil(max_distance / bin_width)
atom_mode = "c3prime" 

def residue_distances(atoms):
    """
    Compute pairwise distances between residues in an RNA structure.

    Parameters
    ----------
    atoms : list
        List of atom entries, where each entry contains chain, residue ID, and 3D coordinates.

    Returns
    -------
    list of tuples
        Each tuple is `(res_i, res_j, distance)` representing the distance between two residues.
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
    Compute counts of residue-residue distances for reference and base pair-specific distributions.

    Parameters
    ----------
    atoms : list
        List of atom entries, where each entry contains chain, residue ID, and 3D coordinates.

    Returns
    -------
    tuple`(reference_counts, pair_counts)` where `reference_counts` is a list of counts for all residues,
    and `pair_counts` is a dictionary of counts per base pair.
    """

    pair_counts = {bp: [0] * num_bins for bp in base_pairs}
    reference_counts = [0] * num_bins

    for res_i, res_j, distance in residue_distances(atoms):

        bin_index = int(distance / bin_width)

        # if bin_index >= num_bins:
        #     continue

        if bin_index >= num_bins:
            bin_index = num_bins - 1

        base_pair = normalize_pair(res_i, res_j)

        reference_counts[bin_index] += 1
        pair_counts[base_pair][bin_index] += 1

    return reference_counts, pair_counts

def frequencies(counts):
    """
    Convert raw counts into normalized frequency values.

    Parameters
    ----------
    counts : list or np.ndarray
        Raw counts per distance bin.

    Returns
    -------
    np.ndarray
        Array of normalized frequencies corresponding to input counts.
    """

    counts = np.asarray(counts, dtype=float)
    total = counts.sum()
    if total == 0:
        return np.zeros_like(counts)
    return counts / total

def score(reference_frequency, pair_frequency):
    """
    Compute the scoring term u_ij = -log(f_ij / f_xx) for a residue pair.

    Parameters
    ----------
    reference_frequency : float
        Reference frequency for a given distance bin.
    pair_frequency : float
        Frequency for the specific residue pair in the same bin.

    Returns
    -------
    float
        Score u_ij; returns `inf` if either frequency is zero.
    """

    if reference_frequency == 0 or pair_frequency == 0:
        return float("inf")
    return -math.log(pair_frequency / reference_frequency)

def residue_distances_all_atom(residues):
    """
    Compute minimal all-atom distances between residue pairs
    within the same chain.

    Parameters
    ----------
    residues : list of tuples - [(chain_id, resseq, resname, coords(N,3)), ...]
        Each element is (chain_id, resseq, resname, coords),
        where coords is a (N_atoms, 3) NumPy array.

    Returns
    -------
    list of tuples
        (resname_i, resname_j, min_distance)
    """

    if not residues:
        return []

    chains = [r[0] for r in residues]
    resnames = [r[2] for r in residues]
    coords_list = [r[3] for r in residues]

    n = len(residues)
    distances = []

    for i in range(n - position_skip):
        chain_i = chains[i]
        coords_i = coords_list[i]

        for j in range(i + position_skip, n):
            if chains[j] != chain_i:
                continue

            coords_j = coords_list[j]

            # min squared distance between two point clouds
            diff = coords_i[:, None, :] - coords_j[None, :, :]
            d2_min = float(np.sum(diff * diff, axis=2).min())

            if d2_min < max_distance_sq:
                distances.append((resnames[i], resnames[j], math.sqrt(d2_min)))

    return distances


def distance_counts_all_atom(residues):
    """
    Compute counts for reference and base pair-specific 
    distributions using all-atom min distances.

    For each residue pair, distances are accumulated into:
    - a global (reference) distribution
    - base-pair–specific distributions

    Parameters
    ----------
    residues : list of tuples
        Each element is (chain_id, resseq, resname, coords),
        where coords is a (N_atoms, 3) NumPy array.

    Returns
    -------
    reference_counts : list of int
        Histogram of all residue-pair distances.
    pair_counts : dict
        Mapping base_pair -> histogram of distances.
    """

    pair_counts = {bp: [0] * num_bins for bp in base_pairs}
    reference_counts = [0] * num_bins

    for res_i, res_j, distance in residue_distances_all_atom(residues):
        bin_index = int(distance / bin_width)
        if bin_index >= num_bins:
            bin_index = num_bins - 1

        bp = normalize_pair(res_i, res_j)
        reference_counts[bin_index] += 1
        pair_counts[bp][bin_index] += 1

    return reference_counts, pair_counts

def save_pairwise_distances(structure_dir="data/structures/test", output_file="data/distances.csv"):
    """
    Extract all base pair distances from RNA structures in a directory
    and save them in a CSV/TSV file.

    Parameters
    ----------
    structure_dir : str
        Folder containing PDB/CIF files.
    output_file : str
        Path to the CSV file to save distances.
    """

    if not os.path.isdir(structure_dir):
        raise FileNotFoundError(f"Directory not found: {structure_dir}")

    structure_files = [
        os.path.join(structure_dir, f)
        for f in os.listdir(structure_dir)
        if f.lower().endswith((".pdb", ".cif", ".mmcif"))
    ]

    if not structure_files:
        raise RuntimeError(f"No PDB/CIF files found in {structure_dir}")

    with open(output_file, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["structure_file", "base_pair", "distance"])

        for struct_file in structure_files:
            if atom_mode == "all_atom":
                residues = extract_all_atom_residues(struct_file)
                distances = residue_distances_all_atom(residues)
            else:
                atoms = extract_c3_atoms(struct_file)
                distances = residue_distances(atoms)

            for res_i, res_j, d in distances:
                bp = normalize_pair(res_i, res_j)
                writer.writerow([os.path.basename(struct_file), bp, d])

    print(f"Distances saved to {output_file}")
