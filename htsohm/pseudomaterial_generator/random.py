from random import random

import numpy as np
from Crypto.Hash import SHA256
from zlib import crc32

from htsohm import config
from htsohm.db import Material
from htsohm.structure import Structure

# random number generator
def bytes_to_float(b):
    return float(crc32(b) & 0xffffffff) / 2**32

def str_to_float(s, encoding="utf-8"):
    return bytes_to_float(s.encode(encoding))

def random_number(n):
    # use SHA256 to get random number from some input value
    h = SHA256.new('{}'.format(n).encode())
    return str_to_float(h.hexdigest())

# uniform selection from random number
def uniform_selection(a, b, r):
    return b - r * (b - a)

def random_number_density(number_density_limits, structure, r): #lattice_constants):
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
    v = structure.volume()
    min_atoms = int(min_ND * v)
    max_atoms = int(max_ND * v)
    if min_atoms < 2:
        min_atoms = int(2)
    return int(uniform_selection(min_atoms, max_atoms, r))

def generate_material(run_id, seed, config):
    """Create records for pseudomaterial simulation and structure data."

    Args:
        run_id (str): identification string for run.
        number_of_atomtypes (int): number of different chemical species used to
            populate the unit cell.

    Returns:
        material (sqlalchemy.orm.query.Query): database row for storing 
            simulation data specific to the material. See
            `htsohm/db/material.py` for more information.

    """
    number_of_atom_types    = config["number_of_atom_types"]
    lattice_limits          = config["lattice_constant_limits"]
    number_density_limits   = config["number_density_limits"]
    epsilon_limits          = config["epsilon_limits"]
    sigma_limits            = config["sigma_limits"]
    max_charge              = config["charge_limit"]

    ########################################################################
    # create database row
    material = Material(run_id)
    material.seed = seed

    # initialize incrementable value for random number generation
    r = material.seed

    # create structure object
    structure = Structure(material.uuid)

    # assign lattice constants
    structure.lattice_constants['a'] = uniform_selection(*lattice_limits, random_number(r))
    r += 1      # increment value to generate a different, reproducible random number
    structure.lattice_constants['b'] = uniform_selection(*lattice_limits, random_number(r))
    r += 1
    structure.lattice_constants['c'] = uniform_selection(*lattice_limits, random_number(r))
    r += 1

    # store unit cell volume to row
    material.ap_unit_cell_volume = structure.volume()
    
    # assign Lennard-Jones parameters
    for chemical_id in range(number_of_atom_types):
        structure.atom_types.append(
                {"chemical-id" : "A_{}".format(chemical_id),
                 "sigma"       : uniform_selection(*sigma_limits, random_number(r)),
                 "epsilon"     : uniform_selection(*epsilon_limits, random_number(r + 1))})
        r += 2
    
    # calculate random number of atom-sites
    number_of_atoms = random_number_density(number_density_limits, structure, random_number(r))
    r += 1

    # store number density to row
    material.ap_number_density = number_of_atoms / structure.volume()

    # assign atom-site positions and calculate avg. sigma/epsilon values
    sigma_sum, epsilon_sum = 0, 0
    for i in range(number_of_atoms):
        chemical_id = int(uniform_selection(0, number_of_atom_types - 1, random_number(r)))
        sigma_sum += structure.atom_types[chemical_id]["sigma"]
        epsilon_sum += structure.atom_types[chemical_id]["epsilon"]
        structure.atom_sites.append(
                {"chemical-id" : "A_{}".format(chemical_id),
                 "x-frac" : random_number(r + 1),
                 "y-frac" : random_number(r + 2),
                 "z-frac" : random_number(r + 3),
                 "charge" : 0.0})
        r += 4
    material.ap_average_sigma = sigma_sum / number_of_atoms
    material.ap_average_epsilon = epsilon_sum / number_of_atoms
    
    # assign atom-site partial charges
    for i in range(len(structure.atom_sites)):
        a0 = max_charge - abs(structure.atom_sites[i]["charge"])
        j = int(uniform_selection(0, len(structure.atom_sites) - 1, random_number(r)))
        a1 = max_charge - abs(structure.atom_sites[j]["charge"])
        dq = uniform_selection(0, min([a0, a1]), random_number(r + 1))
        structure.atom_sites[i]["charge"] += dq
        structure.atom_sites[j]["charge"] -= dq
        r += 2

    return material, structure

def new_material(run_id, config):
    seed = random()
    return generate_material(run_id, seed, config)
