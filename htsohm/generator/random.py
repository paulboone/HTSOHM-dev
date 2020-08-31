from math import ceil
from random import choice, random, randrange, uniform, randint

from htsohm.db import Material, Structure, AtomTypes, AtomSite
from htsohm.max_pair_distance import min_pair_distance

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

def new_material(config, attempts=10):
    for _ in range(attempts):
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

        if structure.min_pair_distance * structure.a > config['minimum_site_distance']:
            return Material(structure=structure, number_density=number_of_atoms / structure.volume)

    raise(Exception("Failed to create a material that satisfied min site distance requirements in allowed number of attempts"))


def find_atom_site_with_minimum_distance(current_positions, distance, uc_a, num_trials=100):
    for _ in range(num_trials):
        trial_pos = (random(), random(), random())
        if min_pair_distance(current_positions + [trial_pos]) * uc_a > distance:
            return trial_pos
    return None
