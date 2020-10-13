from itertools import product

import numpy as np
from numpy.random import choice

def choose_parents(num_parents, ids, props, bin_materials, r=1):
    all_parent_bins = [coords for coords, mats in np.ndenumerate(bin_materials) if len(mats) > 0]

    imax = len(bin_materials)
    jmax = len(bin_materials[0])

    eligible_parent_bins = []
    for coords in all_parent_bins:
        ci, cj = coords
        ilo = max(ci - r, 0)
        ihi = min(ci + r, imax - 1)
        jlo = max(cj - r, 0)
        jhi = min(cj + r, jmax - 1)
        neighbor_count = 0
        for i,j in product(range(ilo, ihi + 1), range(jlo, jhi + 1)):
            if len(bin_materials[i][j]) == 0:
                eligible_parent_bins.append(coords)
                break

    parent_bins = choice(np.array(eligible_parent_bins, "i,i"), num_parents)
    parent_indices = [choice(bin_materials[bin[0]][bin[1]], 1)[0] for bin in parent_bins]

    return [ids[i] for i in parent_indices], [props[i] for i in parent_indices]
