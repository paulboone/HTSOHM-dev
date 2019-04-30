import numpy as np

from htsohm.db import Material
from numpy.random import choice

def choose_parents(num_parents, box_d, box_range, bin_materials, score_by_empty_neighbors=False):
    raw_bins = [(i, len(mats)) for i, mats in np.ndenumerate(bin_materials) if len(mats) > 0]
    if not score_by_empty_neighbors:
        bins = raw_bins
    else:
        bins = []
        r = 2
        imax = len(bin_materials)
        jmax = len(bin_materials[0])
        for (ci, cj), mat_count in raw_bins:
            ilo = max(ci - r, 0)
            ihi = min(ci + r, imax - 1)
            jlo = max(cj - r, 0)
            jhi = min(cj + r, jmax - 1)
            neighbor_count = 0
            for i in range(ilo, ihi + 1):
                for j in range(jlo, jhi + 1):
                    if len(bin_materials[i][j]) == 0:
                        neighbor_count += 1
            bins.append([(ci, cj), 25 - neighbor_count])

    bins.sort(key=lambda x: x[1])

    # since we are sorted, this is the bin with the most materials
    cutoff_index = num_parents - 1 if num_parents - 1 < len(bins) else -1
    cutoff = bins[cutoff_index][1]

    # limit to ALL materials that are within the cutoff. This is necessary because our weighting is
    # based on an integer value here, as opposed to the float values for the convex_hull methods.
    bins = np.array([x for x in bins if x[1] <= cutoff])

    bin_indices = bins[:, 0]
    bin_weights = bins[:, 1].astype(float)

    # calculate weights by subtracting the # materials per bin from the total weight to get a
    bin_weights = bin_weights.sum() / bin_weights
    bin_weights /= bin_weights.sum()
    parent_bins = choice(bin_indices, num_parents, p=bin_weights)

    parent_indices = []
    for bin in parent_bins:
        parent_indices.append(choice(bin_materials[bin[0]][bin[1]], 1)[0])

    return [box_d[i] for i in parent_indices], [box_range[i] for i in parent_indices], bins
