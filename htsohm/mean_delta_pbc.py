
import numpy as np

def minimum_distance_v(v1, v2):
    vdiff = v2 - v1
    if abs(vdiff) < 0.5:
        return vdiff
    elif v1 > 0.5:
        return 1 - abs(vdiff)
    else:
        return -1 + abs(vdiff)


def minimum_distance_point(p1, p2):
    return [minimum_distance_v(v1, v2) for (v1, v2) in zip(p1, p2)]

def center_delta_pbc(points):
    return mean_delta_pbc(points)

def mean_delta_pbc(points):
    """ calculates the best mean for a group of points through periodic boundary conditions with the
        bounds of {0,1} for all axes. It does this by calculating the minimum offset from the first
        point to all the other points, either within the unit cell, or through any boundaries. Then
        we calculate the mean relative to the first point to get the mean of the cluster, which can
        then be added to the absolute coordinates of the first point to get the mean in absolute
        point space (after we modulus by 1.)

        Warning: this may not be mathematically robust! But the point is for this to be just a rough
        metric of how close a pseudomaterial is to having all the points overlap (the mean delta of
        this case would be 0).
    """
    means = []
    distances = []
    deltas = []

    best_range = None
    best_mean = None
    for p0 in points:
        ps = [minimum_distance_point(p0, p1) for p1 in points]
        rel_mean = np.mean(ps, axis=0)
        deltas = np.array(ps) - rel_mean
        min_range = max(np.sum(deltas ** 2, axis=1) ** 0.5)
        mean = (rel_mean + p0) % 1.0

        if (best_range is None) or (min_range < best_range):
            best_range = min_range
            best_mean = mean

    return best_mean, best_range
