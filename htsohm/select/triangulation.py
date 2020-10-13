
from random import random

import numpy as np
from numpy.random import choice
from scipy.spatial import Delaunay

def choose_parents_hull(triang, props, num_parents):
    hull_point_indices = np.unique(triang.convex_hull.flatten())
    point_weights = {i:0.0 for i in hull_point_indices}
    
    for edge in triang.convex_hull:
        distance = np.sqrt(np.sum((props[edge[0]] - props[edge[1]]) ** 2))
        point_weights[edge[0]] += distance
        point_weights[edge[1]] += distance

    point_weight_arr = [[point_weights[i], i] for i in point_weights.keys()]
    point_weight_arr.sort(key=lambda x: x[0])
    point_weight_arr = np.array(point_weight_arr)[-num_parents:]

    total_weight = point_weight_arr[:,0].sum()
    point_weight_arr[:, 0] /= total_weight

    parent_indices = choice(point_weight_arr[:,1], num_parents, p=point_weight_arr[:,0])
    parent_indices.sort()
    parent_indices = [int(i) for i in parent_indices]

    return parent_indices

def triangle_area(p1, p2, p3):
    return (p1[0]*(p2[1] - p3[1]) + p2[0]*(p3[1] - p1[1]) + p3[0]*(p1[1] - p2[1])) / 2

def choose_parents_simplices(triang, props, num_parents, num_best_triangles):
    areas = [[triangle_area(props[p1], props[p2], props[p3]), p1, p2, p3] for p1, p2, p3 in triang.simplices]
    num_best_triangles = min(len(areas), num_best_triangles) # necessary for small generations
    areas.sort()
    areas = np.array(areas[-num_best_triangles:])
    total_area = areas[:,0].sum()
    areas[:,0] /= total_area
    triangle_indices = choice(range(num_best_triangles), num_parents, p=areas[:,0])
    parent_indices = [int(areas[t, choice([0,1,2]) + 1]) for t in triangle_indices]
    return parent_indices


def choose_parents(num_parents, ids, props, simplices_or_hull):
    triang = Delaunay(props)
    if simplices_or_hull == 'simplices':
        parent_indices = choose_parents_simplices(triang, props, num_parents, num_parents)
    elif simplices_or_hull == 'hull':
        parent_indices = choose_parents_hull(triang, props, num_parents)
    else:
        raise(Exception("simplices_or_hull must be defined as 'simplices' or 'hull'"))

    return [ids[i] for i in parent_indices], [props[i] for i in parent_indices]
