from math import ceil
from random import choice, random, randrange, uniform, randint

from htsohm.db import Material, Structure, LennardJones, AtomSite

def random_atom_types(num_atom_types, config):
    return [LennardJones(
        sigma = uniform(*config["sigma_limits"]),
        epsilon = uniform(*config["epsilon_limits"])
    ) for i in range(num_atom_types)]

def random_atom_sites(num_sites, atom_types):
    return [AtomSite(
        lennard_jones = atom_types[choice(range(len(atom_types)))],
        x = random(), y = random(), z = random(),
        q = 0.0
    ) for i in range(num_sites)]

def new_material(config):
    structure=Structure(
        a = uniform(*config["lattice_constant_limits"]),
        b = uniform(*config["lattice_constant_limits"]),
        c = uniform(*config["lattice_constant_limits"]),
        lennard_jones = random_atom_types(config["number_of_atom_types"], config),
    )

    if "fix_atoms" in config:
        number_of_atoms = config['fix_atoms']
    else:
        number_of_atoms = random_number_density(config["number_density_limits"], structure.volume)

    number_of_atoms = randint(*config["num_atoms_limits"])
    structure.atom_sites = random_atom_sites(number_of_atoms, structure.lennard_jones)

    return Material(structure=structure, number_density=number_of_atoms / structure.volume)

def random_number_density(number_density_limits, volume):
    """Produces random number for atom-sites in a unit cell,
    constrained by some number density limits and a volume"""

    min_atoms = ceil(number_density_limits[0] * volume)
    max_atoms = int(number_density_limits[1] * volume)
    if min_atoms < 1:
        min_atoms = 1
    if min_atoms >= max_atoms:
        return min_atoms

    return randrange(min_atoms, max_atoms + 1, 1)
