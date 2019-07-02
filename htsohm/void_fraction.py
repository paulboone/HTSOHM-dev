
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
    lattice_indices = itertools.product(range(xi_max), range(yi_max), range(zi_max))
    lattice_fill = np.zeros((xi_max, yi_max, zi_max))

    for x, y, z, r in atoms:
        for lxi, lyi, lzi in lattice_indices:
            lx = lxi * dx; ly = lyi * dy; lz = lzi * dz
            if math.sqrt((lx - x)**2 + (ly - y)**2 + (lz -z)**2) < r:
                lattice_fill[lxi, lyi, lzi] = 1.0

    return 1 - np.sum(lattice_fill) / (xi_max * yi_max * zi_max)

def get_radial_distances(xyz_values, cutoff, calc_product=True):
    # gets list of radial distances to all framework atoms in one direction of 3d space
    if calc_product:
        all_atoms = itertools.product(xyz_values, repeat=3)
    else:
        all_atoms = xyz_values

    all_atoms_r = [(a,np.sqrt(np.sum(np.array(a)**2))) for a in all_atoms]
    all_atoms_r.sort(key=lambda x: x[1])
    all_atoms_r = list(filter(lambda x: x[1]<cutoff, all_atoms_r))
    return [x[1] for x in all_atoms_r]
