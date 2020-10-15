import numpy as np

from htsohm.db import Material
from numpy.random import choice


def choose_parent_bins_from_weighted_bin_list(bins, num_parents):
    """
    bins: (bin_tuple, weight)

    The weights of all the bins will be summed up and each bin will be normalized.
    """

    bins.sort(key=lambda x: x[1])

    # since we are sorted, this is the bin with the most materials
    cutoff_index = num_parents - 1 if num_parents - 1 < len(bins) else -1
    cutoff = bins[cutoff_index][1]

    # limit to ALL materials that are within the cutoff. This is necessary because our weighting is
    # based on an integer value here, as opposed to the float values for the convex_hull methods.
    bins = np.array([x for x in bins if x[1] <= cutoff], dtype=object)

    bin_tuples  = bins[:, 0]
    bin_weights = bins[:, 1].astype(float)

    # calculate weights by subtracting the # materials per bin from the total weight to get a
    bin_weights = bin_weights.sum() / bin_weights
    bin_weights /= bin_weights.sum()

    return choice(bin_tuples, num_parents, p=bin_weights)


def choose_parents(num_parents, ids, props, bin_materials):
    bins = [(bin_tuple, len(mats)) for bin_tuple, mats in np.ndenumerate(bin_materials) if len(mats) > 0]
    parent_bin_tuples = choose_parent_bins_from_weighted_bin_list(bins, num_parents)
    parent_indices = [choice(bin_materials[bin_tuple], 1)[0] for bin_tuple in parent_bin_tuples]
    return [ids[i] for i in parent_indices], [props[i] for i in parent_indices]
