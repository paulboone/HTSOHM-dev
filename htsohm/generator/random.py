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
        lennard_jones = choice(atom_types),
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

    number_of_atoms = randint(*config["num_atoms_limits"])
    structure.atom_sites = random_atom_sites(number_of_atoms, structure.lennard_jones)

    return Material(structure=structure, number_density=number_of_atoms / structure.volume)
