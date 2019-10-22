import numpy as np
from numpy.random import choice

from htsohm.select.density_bin import choose_parent_bins_from_weighted_bin_list
from htsohm.db import Material


def choose_parents(num_parents, box_d, box_range, bin_materials):
    raw_bins = [(i, len(mats)) for i, mats in np.ndenumerate(bin_materials) if len(mats) > 0]
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

    parent_bins = choose_parent_bins_from_weighted_bin_list(bins, num_parents)
    parent_indices = [choice(bin_materials[bin[0]][bin[1]], 1)[0] for bin in parent_bins]

    return [box_d[i] for i in parent_indices], [box_range[i] for i in parent_indices]
