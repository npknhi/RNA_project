import math
from model.pair import set_pairs, normalize_pair 

# Parameters
nucleotides = ["A", "U", "G", "C"]
base_pairs = set_pairs(nucleotides)
max_distance = 20
position_skip = 4
maximum_score = 10


def euclidian_distance(x1, y1, z1, x2, y2, z2):
    """
    Compute Euclidean distance between two 3D points.
    """
    return math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2 + (z2 - z1) ** 2)


def residue_distances(atoms):
    """
    Compute distances between residues within the same chain,
    skipping positions according to position_skip.
    
    Parameters:
        atoms (list of tuples): Each tuple is (chain_id, residue_name, x, y, z)
    
    Returns:
        list of tuples: (residue_i, residue_j, distance)
    """
    distances = []
    for i in range(len(atoms) - position_skip):
        chain_i, residue_i, x_i, y_i, z_i = atoms[i]

        for j in range(i + position_skip, len(atoms)):
            chain_j, residue_j, x_j, y_j, z_j = atoms[j]
            if chain_j != chain_i:
                continue
            distance = euclidian_distance(x_i, y_i, z_i, x_j, y_j, z_j)
            if distance < max_distance:
                distances.append((residue_i, residue_j, distance))
    return distances


def distance_counts(atoms):
    """
    Compute counts of atomic distances.
    
    Returns:
        tuple: (reference_counts, pair_counts)
        - reference_counts: list of counts for all pairs
        - pair_counts: dict mapping base_pair -> list of counts
    """
    # Initialize counts
    pair_counts = {bp: [0] * max_distance for bp in base_pairs}
    reference_counts = [0] * max_distance

    # Count each distance
    for residue_i, residue_j, distance in residue_distances(atoms):
        distance_bin = int(math.floor(distance))
        base_pair = normalize_pair(residue_i, residue_j)
        reference_counts[distance_bin] += 1
        pair_counts[base_pair][distance_bin] += 1

    return reference_counts, pair_counts


def frequencies(counts):
    """
    Normalize counts to frequencies.
    """
    total = sum(counts)
    if total == 0:
        return [0] * len(counts)
    return [c / total for c in counts]


def score(reference_frequency, pair_frequency):
    """
    Compute the score of a pair of nucleotides:
    u_ij = -log(f_ij / f_xx)
    """
    if reference_frequency == 0 or pair_frequency == 0:
        return float("inf")
    return -math.log(pair_frequency / reference_frequency)
