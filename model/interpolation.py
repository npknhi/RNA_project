def linear_interpolate(profile, d, bin_size=1.0):
    if d >= len(profile):
        return profile[-1]

    i = int(d // bin_size)
    t = (d % bin_size) / bin_size

    if i >= len(profile) - 1:
        return profile[-1]

    return profile[i] * (1 - t) + profile[i+1] * t
