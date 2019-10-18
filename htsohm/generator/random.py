from random import choice, random, randrange, uniform

from htsohm.db import Material, Structure, LennardJones, AtomSite
from htsohm.generator.utilities import random_number_density

def new_material(config):
    # randomly generate atom types
    lj_atom_types = [LennardJones(
        sigma = uniform(*config["sigma_limits"]),
        epsilon = uniform(*config["epsilon_limits"])
    ) for i in range(config["number_of_atom_types"])]

    # randomly generate structure
    structure=Structure(
        a = uniform(*config["lattice_constant_limits"]),
        b = uniform(*config["lattice_constant_limits"]),
        c = uniform(*config["lattice_constant_limits"]),
        lennard_jones = lj_atom_types
    )

    # assign atom-site positions
    if "fix_atoms" in config:
        number_of_atoms = config['fix_atoms']
    else:
        number_of_atoms = random_number_density(config["number_density_limits"], structure)
    
    for i in range(number_of_atoms):
        atom_type_index = choice(range(config["number_of_atom_types"]))
        structure.atom_sites.append(
                AtomSite(
                    lennard_jones = lj_atom_types[atom_type_index],
                    x = random(), y = random(), z = random(),
                    q = 0.0))

    m = Material(structure=structure, number_density=number_of_atoms / structure.volume)

    return m
