def linear_interpolation(x0, y0, x1, y1, x):
    """
    Perform linear interpolation between two known points (x0, y0) and (x1, y1).
    
    ---
    Arguments:
        x0 (float): x-coordinate of the first point
        y0 (float): y-coordinate of the first point
        x1 (float): x-coordinate of the second point
        y1 (float): y-coordinate of the second point
        x (float): x-coordinate where we want to estimate the value
    ---
    Returns:
        float: interpolated y value at position x
    """
    
    if x1 == x0:
        return y0
    return y0 + (y1 - y0) * (x - x0) / (x1 - x0)