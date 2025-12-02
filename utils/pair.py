from itertools import combinations_with_replacement

def normalize_pair(b1, b2):
    return "".join(sorted([b1, b2]))

def set_pairs(nucleotides):
    """
    Return all unique sorted nucleotide pairs (including self-pairs)
    """
    nucleotides = sorted(nucleotides)
    return ["".join(pair) for pair in combinations_with_replacement(nucleotides, 2)]
