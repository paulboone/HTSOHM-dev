
import itertools
import math

import numpy as np

def calculate_void_fraction(atoms, box, points_per_angstrom=10):
    """ assumes box starts at (0,0,0)
    """
    xi_max = box[0] * points_per_angstrom
    yi_max = box[1] * points_per_angstrom
    zi_max = box[2] * points_per_angstrom

    dx = 1 / points_per_angstrom
    dy = 1 / points_per_angstrom
    dz = 1 / points_per_angstrom

    lattice_fill = np.zeros((xi_max, yi_max, zi_max))

    for x, y, z, r in atoms:
        lattice_indices = itertools.product(range(xi_max), range(yi_max), range(zi_max))
        for lxi, lyi, lzi in lattice_indices:
            lx = lxi * dx; ly = lyi * dy; lz = lzi * dz
            if math.sqrt((lx - x)**2 + (ly - y)**2 + (lz -z)**2) < r:
                lattice_fill[lxi, lyi, lzi] = 1.0

    return 1 - np.sum(lattice_fill) / (xi_max * yi_max * zi_max)
