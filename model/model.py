import numpy as np
from .pair import normalize_pair
from .interpolation import linear_interpolate

class RNAModel:
    def __init__(self, profile_folder):
        self.profiles = {}
        self.load_profiles(profile_folder)

    def load_profiles(self, folder):
        pairs = ["AA","AU","AC","AG","UU","UC","UG","CC","CG","GG"]
        for p in pairs:
            path = f"{folder}/{p}.txt"
            self.profiles[p] = np.loadtxt(path)

    def score_distance(self, b1, b2, d):
        pair = normalize_pair(b1, b2)
        if pair not in self.profiles:
            return 0.0
        return linear_interpolate(self.profiles[pair], d)
