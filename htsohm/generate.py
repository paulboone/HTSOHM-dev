# standard library imports
import os
from math import floor
from random import choice, random, randrange, uniform

# it stat third party imports
from functools import reduce
import numpy as np
import yaml

# local application/library specific imports
from htsohm.runDB_declarative import session, RunData
import htsohm.utilities as utl

def write_material_config(run_id):
    """ Write material-parameters to run-configuration file.
    
    The parameters written by this function define the limits for different values written to the
    structure and forcefield definition files for RASPA. Among the limits defined here are crystal
    lattice constants, number density, partial atomic charges, and Lennard-Jones parameters (sigma
    and epsilon).
    """
    wd = os.environ['HTSOHM_DIR']      # specify $HTSOHM_DIR as working directory
    config_file = os.path.join(wd, 'config', run_id + '.yaml')
    with open(config_file) as file:
        run_config = yaml.load(file)

    run_config.update({
        "number-density-limits"     : [0.000013907, 0.084086],
        "lattice-constant-limits"   : [13.098, 52.392],
        "epsilon-limits"            : [1.258, 513.264],
        "sigma-limits"              : [1.052, 6.549],
        "charge-limit"              : 0.,
        "elemental-charge"          : 0.0001
    })

    with open(config_file, "w") as file:
        yaml.dump(run_config, file, default_flow_style=False)

    return run_config

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

def write_seed_definition_files(run_id, number_of_materials, number_of_atomtypes):
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

    material_config         = write_material_config(run_id)
    lattice_limits          = material_config["lattice-constant-limits"]
    number_density_limits   = material_config["number-density-limits"]
    epsilon_limits          = material_config["epsilon-limits"]
    sigma_limits            = material_config["sigma-limits"]
    max_charge              = material_config["charge-limit"]
    elem_charge             = material_config["elemental-charge"]

    wd = os.environ['HTSOHM_DIR']      # specify $HTSOHM_DIR as working directory
    ff_dir = os.environ['FF_DIR']      # output force-field files to $FF_DIR
    mat_dir = os.environ['MAT_DIR']    # output .cif-files to $MAT_DIR

    materials = session.query(RunData).filter(RunData.run_id == run_id, RunData.generation == 0).all()
    for material in materials:           # each iteration creates a new material
        material_name = run_id + '-' + str(material.id)        # this will be replaced with primary_key

        def_dir = os.path.join(ff_dir, material_name)       # directory for material's force field
        os.mkdir(def_dir)
        force_field_file = os.path.join(def_dir, 'force_field.def')      # for overwriting LJ-params
        utl.write_force_field(force_field_file)

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
        utl.write_mixing_rules(mix_file, atom_types)
        psu_file = os.path.join(def_dir, 'pseudo_atoms.def')             # define atom-types
        utl.write_pseudo_atoms(psu_file, atom_types)

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
        utl.write_cif_file(cif_file, lattice_constants, atom_sites)
