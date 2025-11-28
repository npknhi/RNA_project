BASES = ["A", "U", "C", "G"]

def normalize_pair(b1, b2):
    return "".join(sorted([b1, b2]))
