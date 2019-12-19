from math import ceil
from random import choice, random, randrange, uniform, randint

from htsohm.db import Material, Structure, AtomTypes, AtomSite

def random_atom_types(num_atom_types, config):
    return [AtomTypes(
        sigma = uniform(*config["sigma_limits"]),
        epsilon = uniform(*config["epsilon_limits"])
    ) for i in range(num_atom_types)]

def random_atom_sites(num_sites, atom_types):
    return [AtomSite(
        atom_types = choice(atom_types),
        x = random(), y = random(), z = random(),
        q = 0.0
    ) for i in range(num_sites)]

def new_material(config):
    a = uniform(*config["lattice_constant_limits"])
    if config["lattice_cubic"]:
        b = a
        c = a
    else:
        b = uniform(*config["lattice_constant_limits"])
        c = uniform(*config["lattice_constant_limits"])

    structure=Structure(
        a = a, b = b, c = c,
        atom_types = random_atom_types(config["number_of_atom_types"], config),
    )

    number_of_atoms = randint(*config["num_atoms_limits"])
    structure.atom_sites = random_atom_sites(number_of_atoms, structure.atom_types)

    return Material(structure=structure, number_density=number_of_atoms / structure.volume)
