import os

from random import choice, random, randrange, uniform
from functools import reduce
from math import fsum
import numpy as np
from math import floor
import yaml

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
    max_number_density = number_density_limits[1]
    a = lattice_constants["a"]
    b = lattice_constants["b"]
    c = lattice_constants["c"]
    max_number_of_atoms = int(max_number_density * a * b * c)
    number_of_atoms = randrange(2, max_number_of_atoms, 1)
    number_density = round(number_of_atoms / (a * b * c))
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
    seed = session.query(RunData).filter(RunData.run_id == run_id, RunData.generation == 0).all()
    seed_ids = []
    for material in seed:              # get list of seed-material IDs
        seed_ids.append(material.id)

    material_config = write_material_config(run_id)
    lattice_limits          = material_config["lattice-constant-limits"]
    number_density_limits   = material_config["number-density-limits"]
    epsilon_limits          = material_config["epsilon-limits"]
    sigma_limits            = material_config["sigma-limits"]
    max_charge              = material_config["charge-limit"]
    elem_charge             = material_config["elemental-charge"]

    wd = os.environ['HTSOHM_DIR']      # specify $HTSOHM_DIR as working directory
    ff_dir = os.environ['FF_DIR']      # output force-field files to $FF_DIR
    mat_dir = os.environ['MAT_DIR']    # output .cif-files to $MAT_DIR

    materials = []
    for material in seed_ids:           # each iteration creates a new material
        material_id = run_id + '-' + str(material)        # this will be replaced with primary_key

        def_dir = os.path.join(ff_dir, material_id)       # directory for material's force field
        os.mkdir(def_dir)
        force_field_file = os.path.join(def_dir, 'force_field.def')      # for overwriting LJ-params
        utl.write_force_field(force_field_file)

        ########################################################################
        # define pseudo atom types by randomly-generating sigma and epsilon values
        atom_types = []
        for chemical_id in range(number_of_atomtypes):
            epsilon   = round(uniform(*epsilon_limits), 4)
            sigma     = round(uniform(*sigma_limits), 4)
            charge    = 0.             # charge assignment to be re-implemented!!!
            atom_type = {
                "chemical-id" : "A_%s" % chemical_id,
                "charge"      : charge,
                "epsilon"     : epsilon,
                "sigma"       : sigma}
            atom_types.append(atom_type)

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
            atom_type = choice(atom_types)
            atom_site = {
                "chemical-id" : atom_type["chemical-id"],
                "x-frac"      : round(random(), 4),
                "y-frac"      : round(random(), 4),
                "z-frac"      : round(random(), 4)}
            atom_sites.append(atom_site)

        cif_file = os.path.join(mat_dir, material_id + ".cif")           # structure file
        utl.write_cif_file(cif_file, lattice_constants, atom_sites)
