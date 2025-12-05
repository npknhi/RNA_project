from itertools import combinations_with_replacement

def normalize_pair(b1, b2):
    """
    Return a normalized representation of a nucleotide pair by sorting the two bases.

    Parameters
    ----------
    b1 : str
        First nucleotide.
    b2 : str
        Second nucleotide.

    Returns
    -------
    str
        A string of the two nucleotides in sorted order.
    """
    
    return "".join(sorted([b1, b2]))

def set_pairs(nucleotides):
    """
    Generate all unique sorted nucleotide pairs, including self-pairs.

    Parameters
    ----------
    nucleotides : list
        List of nucleotide characters.

    Returns
    -------
    list of str
        List of strings representing all unique nucleotide pairs.
    """

    nucleotides = sorted(nucleotides)
    return ["".join(pair) for pair in combinations_with_replacement(nucleotides, 2)]
