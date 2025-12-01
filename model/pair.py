from itertools import combinations_with_replacement

def normalize_pair(b1, b2):
    return "".join(sorted([b1, b2]))

def set_pairs(nucleotides):
    nucleotides = sorted(nucleotides)
    return ["".join(pair) for pair in combinations_with_replacement(nucleotides, 2)]
