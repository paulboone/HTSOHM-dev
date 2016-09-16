# standard library imports
from functools import reduce
from math import floor
import os
from random import choice, random, randrange, uniform
from uuid import uuid4

# related third party imports
import numpy as np
import yaml

# local application/library specific imports
from htsohm.runDB_declarative import session, Material
from htsohm.utilities import update_config_file, write_force_field, write_cif_file
from htsohm.utilities import write_mixing_rules, write_pseudo_atoms

def random_number_density(number_density_limits, lattice_constants):
    """Returns some random number of atoms per unit cell, within the defined density limits.
    There are two inputs :
        - number_density limits        a list of the form [minimum, maximum]
        - lattice_constants            a dictionary of crystal lattice parameters of the form
                                           {"a" : <float>,
                                            "b" : <float>,
                                            "c" : <float>}
    And one output:
        - number_of_atoms              some random number of atoms per unit cell, according to
                                       predefined limits.
    """
    max_number_density = number_density_limits[1]
    a = lattice_constants["a"]
    b = lattice_constants["b"]
    c = lattice_constants["c"]
    max_number_of_atoms = int(max_number_density * a * b * c)
    number_of_atoms = randrange(2, max_number_of_atoms, 1)
    return number_of_atoms

def write_seed_definition_files(run_id, number_of_atomtypes):
    """Write .def and .cif files for a randomly-generated porous material.

    Each material is defined by it's structural information (stored in a .cif-file) and force field
    definition files:
    - <material_name>.cif            contains structural information including crystal lattice
                                     parameters and atom-site positions (and corresponding chemical
                                     species).
    - force_field.def                this file can be used to overwrite previously-defined
                                     interactions. by default there are no exceptions.
    - force_field_mixing_rules.def   this file contains sigma and epsilon values to define Lennard-
                                     Jones type interactions.
    - pseudo_atoms.def               this file contains pseudo-atom definitions, including partial
                                     charge, atomic mass, atomic radii, and more.
    """

    material_config         = update_config_file(run_id)
    lattice_limits          = material_config["lattice-constant-limits"]
    number_density_limits   = material_config["number-density-limits"]
    epsilon_limits          = material_config["epsilon-limits"]
    sigma_limits            = material_config["sigma-limits"]
    max_charge              = material_config["charge-limit"]
    elem_charge             = material_config["elemental-charge"]

    wd = os.environ['HTSOHM_DIR']                    # specify $HTSOHM_DIR as working directory
    ff_dir = os.environ['FF_DIR']                    # output force-field files to $FF_DIR
    mat_dir = os.environ['MAT_DIR']                  # output .cif-files to $MAT_DIR

    ########################################################################
    material = Material(run_id, 'none')
    material.seed = True
    material_name = run_id + '-' + material.uuid

    def_dir = os.path.join(ff_dir, material_name)       # directory for material's force field
    os.mkdir(def_dir)
    force_field_file = os.path.join(def_dir, 'force_field.def')      # for overwriting LJ-params
    write_force_field(force_field_file)

    ########################################################################
    # define pseudo atom types by randomly-generating sigma and epsilon values
    atom_types = []
    for chemical_id in range(number_of_atomtypes):
        atom_types.append({
            "chemical-id" : "A_%s" % chemical_id,
            "charge"      : 0.,    # charge assignment to be re-implemented!!!,
            "epsilon"     : round(uniform(*epsilon_limits), 4),
            "sigma"       : round(uniform(*sigma_limits), 4)
        })

    mix_file = os.path.join(def_dir, 'force_field_mixing_rules.def') # LJ-parameters
    write_mixing_rules(mix_file, atom_types)
    psu_file = os.path.join(def_dir, 'pseudo_atoms.def')             # define atom-types
    write_pseudo_atoms(psu_file, atom_types)

    ########################################################################
    # randomly-assign dimensions (crystal lattice constants) and number of atoms per unit cell
    lattice_constants = {"a" : round(uniform(*lattice_limits), 4),
                         "b" : round(uniform(*lattice_limits), 4),
                         "c" : round(uniform(*lattice_limits), 4)}
    number_of_atoms   = random_number_density(number_density_limits, lattice_constants)

    ########################################################################
    # populate unit cell with randomly-positioned atoms of a randomly-selected species
    atom_sites = []
    for atom in range(number_of_atoms):
        atom_sites.append({
            "chemical-id" : choice(atom_types)["chemical-id"],
            "x-frac"      : round(random(), 4),
            "y-frac"      : round(random(), 4),
            "z-frac"      : round(random(), 4)
        })

    cif_file = os.path.join(mat_dir, material_name + ".cif")           # structure file
    write_cif_file(cif_file, lattice_constants, atom_sites)

    material.write_check = 'done'
    return material
