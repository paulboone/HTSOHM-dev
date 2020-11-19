
import itertools
from math import sqrt, ceil

import numpy as np

def calculate_void_fraction(atoms, box, points_per_angstrom=10, probe_r=0.0):
    """calculates a geometric void fraction, given the atom coordinates and diameter and the box size.

    Discretizes the box dimensions into individual cubes where there are points_per_angstrom
    points per angstrom of the box. For each atom in atoms, all cubes within a distance equal to
    the probe radius probe_r and the atom radius (atom diameter / 2) are marked as filled. The
    void fraction is simply the ratio of filled cubes to total cubes.

    Args:
        atoms: an array of tuples (x, y, z, d) containing the atom coordinates x, y, z and the
            atom diameter d.
        box: a tuple containing the length of each side of the box. The box is assumed to start
            at (0,0,0).
        points_per_angstrom: the number of points per angstrom of the box, which determines the
            discretization of the box. The larger the number, the longer the calculations will take to
            perform, and the more accurate the result will be. Default: 10.
        probe_r: the radius of the probe. Default: 0.0 angstroms.

    Returns:
        float: the calculated void fraction.
    """

    xi_max = ceil(box[0] * points_per_angstrom)
    yi_max = ceil(box[1] * points_per_angstrom)
    zi_max = ceil(box[2] * points_per_angstrom)

    dx = box[0] / xi_max
    dy = box[1] / yi_max
    dz = box[2] / zi_max

    lattice_fill = np.zeros((xi_max, yi_max, zi_max))

    for x, y, z, d in atoms:
        r = d/2 + probe_r
        delt_i = ceil(r / dx)
        xi = int(x // dx); yi = int(y // dy); zi = int(z // dz)

        # limit the lattice_indices to only crossing a boundary once in each direction
        lattice_indices = itertools.product(
            range(max(xi - delt_i, -xi_max), min(xi + delt_i + 1, 2*xi_max)),
            range(max(yi - delt_i, -yi_max), min(yi + delt_i + 1, 2*yi_max)),
            range(max(zi - delt_i, -zi_max), min(zi + delt_i + 1, 2*zi_max)))

        for lxi, lyi, lzi in lattice_indices:
            lx = lxi * dx; ly = lyi * dy; lz = lzi * dz
            if sqrt((lx - x)**2 + (ly - y)**2 + (lz -z)**2) < r:
                lattice_fill[lxi % xi_max, lyi % yi_max, lzi % zi_max] = 1.0

    return 1 - np.sum(lattice_fill) / (xi_max * yi_max * zi_max)
