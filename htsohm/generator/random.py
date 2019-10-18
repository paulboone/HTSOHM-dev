from random import choice, random, randrange, uniform

from htsohm.db import Material, Structure, LennardJones, AtomSite
from htsohm.generator.utilities import random_number_density

def new_material(config):
    number_of_atom_types    = config["number_of_atom_types"]
    lattice_limits          = config["lattice_constant_limits"]
    number_density_limits   = config["number_density_limits"]
    epsilon_limits          = config["epsilon_limits"]
    sigma_limits            = config["sigma_limits"]
    max_charge              = config["charge_limit"]

    ########################################################################
    # create database row
    material = Material()
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
    if "fix_atoms" in config:
        number_of_atoms = config['fix_atoms']
    else:
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
                AtomSite(
                    atom_type = "A_{}".format(atom_type),
                    lennard_jones = structure.get_lennard_jones("A_{}".format(atom_type)),
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
