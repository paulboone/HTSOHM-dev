from math import ceil
from random import randrange

def random_number_density(number_density_limits, structure):
    """Produces random number for atom-sites in a unit cell, constrained by
    some number density limits.
    Args:
        number_density_limits (list): min and max number densities, for example:
            [min_(float), max_(float)]
        lattice_constants (dict): crystal lattice constants, for example:
            {"a" : (float),
             "b" : (float),
             "c" : (float)}

    Returns:
        atoms (int): some random number of atom-sites under the imposed limits.
        If the minimum number density results in a unit cell with less than 2
        atom-sites with the given lattice constants, a minimum number density
        of TWO ATOM SITES PER UNIT CELL is imposed.

    """
    min_ND = number_density_limits[0]
    max_ND = number_density_limits[1]
    v = structure.volume
    min_atoms = ceil(min_ND * v)
    max_atoms = int(max_ND * v)
    if min_atoms < 1:
        min_atoms = 1
    if min_atoms >= max_atoms:
        return min_atoms

    return randrange(min_atoms, max_atoms + 1, 1)
