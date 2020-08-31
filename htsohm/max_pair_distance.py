import itertools

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

def max_pair_distance(points):
    """
    does NOT support non-cubic unit cells!
"""
    if len(points) == 1:
        return 0.0
    pairs = np.array([minimum_distance_point(p1, p2) for p1, p2 in itertools.combinations(points, 2)])
    distances = (pairs ** 2).sum(axis=1) ** 0.5
    return max(distances)

def min_pair_distance(points):
    """
    does NOT support non-cubic unit cells!
"""
    if len(points) == 1:
        return 1.0
    pairs = np.array([minimum_distance_point(p1, p2) for p1, p2 in itertools.combinations(points, 2)])
    distances = (pairs ** 2).sum(axis=1) ** 0.5
    return min(distances)
