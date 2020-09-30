from random import choice, random, uniform, randint

from htsohm.db import Material, AtomTypes, AtomSite
from htsohm.pair_distance import min_pair_distance

def random_atom_types(num_atom_types, config):
    return [AtomTypes(
        sigma = uniform(*config["sigma_limits"]),
        epsilon = uniform(*config["epsilon_limits"])
    ) for i in range(num_atom_types)]

def random_charges(num_charges, abs_max_charge):
    """ returns an array containing num_charges charges, such that all charges
    obey the restrictions that abs(charge) <= abs_max_charge and that the total
    charge equals zero.
    """
    total_charge = 0.0
    qs = [0.0] * num_charges
    for i in range(num_charges):
        remaining_assignable_charge = (num_charges - i - 1) * abs_max_charge
        mx = remaining_assignable_charge - total_charge
        mn = -remaining_assignable_charge - total_charge
        qs[i] = uniform(max(mn, -abs_max_charge), min(mx, abs_max_charge))
        total_charge += qs[i]
    return qs

def random_atom_sites(num_sites, atom_types, q=0.0):
    if q == 0.0:
        q = [0.0] * num_sites
    return [AtomSite(
        atom_types = choice(atom_types),
        x = random(), y = random(), z = random(),
        q = q[i]
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


        material=Material(
            a = a, b = b, c = c,
            atom_types = random_atom_types(config["number_of_atom_types"], config),
        )

        number_of_atoms = randint(*config["num_atoms_limits"])

        q = random_charges(number_of_atoms, config["charge_limits"][1])
        material.atom_sites = random_atom_sites(number_of_atoms, material.atom_types, q=q)

        if material.min_pair_distance * material.a > config['minimum_site_distance']:
            return material

    raise(Exception("Failed to create a material that satisfied min site distance requirements in allowed number of attempts"))


def find_atom_site_with_minimum_distance(current_positions, distance, uc_a, num_trials=100):
    for _ in range(num_trials):
        trial_pos = (random(), random(), random())
        if min_pair_distance(current_positions + [trial_pos]) * uc_a > distance:
            return trial_pos
    return None
