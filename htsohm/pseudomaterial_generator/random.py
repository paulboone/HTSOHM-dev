from math import ceil
from random import choice, random, randrange, uniform

from htsohm.db import Material, Structure, LennardJones, AtomSites

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
    min_atoms = ceil(min_ND * v)
    max_atoms = int(max_ND * v)
    if min_atoms < 2:
        min_atoms = 2
    return int(uniform_selection(min_atoms, max_atoms, r))

def generate_material(run_id, config):
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
    material = Material(run_id, None)
    structure = material.structure

    # assign lattice constants
    structure.a = uniform(*lattice_limits)
    structure.b = uniform(*lattice_limits)
    structure.c = uniform(*lattice_limits)

    # store unit cell volume to row
    material.unit_cell_volume = structure.volume
    
    # assign Lennard-Jones parameters
    for i in range(number_of_atom_types):
        structure.lennard_jones.append(
                LennardJones(
                    atom_type  = "A_{}".format(i),
                    sigma      = uniform(*sigma_limits),
                    epsilon    = uniform(*epsilon_limits)))

    # calculate random number of atom-sites
    number_of_atoms = random_number_density(number_density_limits, structure)

    # store number density to row
    material.number_density = number_of_atoms / material.unit_cell_volume

    # assign atom-site positions and calculate avg. sigma/epsilon values
    sigma_sum, epsilon_sum = 0, 0
    for i in range(number_of_atoms):
        # select chemical species
        atom_type = choice(range(number_of_atom_types))

        # sum LJ parameters for averaging later
        sigma_sum += structure.lennard_jones[atom_type].sigma
        epsilon_sum += structure.lennard_jones[atom_type].epsilon

        # set position and add atom-site
        structure.atom_sites.append(
                AtomSites(
                    atom_type = "A_{}".format(atom_type),
                    x = random(), y = random(), z = random(),
                    q = 0.))

    # store avg. sigma/epsilon values to row
    material.average_sigma = sigma_sum / number_of_atoms
    material.average_epsilon = epsilon_sum / number_of_atoms

    # assign atom-site partial charges
    for i in range(number_of_atoms):
        a0 = max_charge - abs(structure.atom_sites[i].q)
        j = choice(range(number_of_atoms))
        a1 = max_charge - abs(structure.atom_sites[j].q)
        dq = float("{0:.6f}".format(uniform(0, min([a0, a1]))))
        structure.atom_sites[i].q += dq
        structure.atom_sites[j].q -= dq

    print("FRAMEWORK NET CHARGE :\t{}".format(sum([e.q for e in structure.atom_sites])))
    return material

def new_material(run_id, config):
    return generate_material(run_id, config)
