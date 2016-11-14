# standard library imports
import sys
import os
from random import choice, random, randrange, uniform
import shutil
from uuid import uuid4

# related third party imports
import numpy as np
import yaml

# local application/library specific imports
from htsohm import config
from htsohm.db import session, Material

def random_number_density(number_density_limits, lattice_constants):
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
    v = lattice_constants["a"] * lattice_constants["b"] * lattice_constants["c"]
    min_atoms = int(min_ND * v)
    max_atoms = int(max_ND * v)
    if min_atoms < 2:
        min_atoms = int(2)
    atoms = randrange(min_atoms, max_atoms, 1)
    return atoms

def write_seed_definition_files(run_id, number_of_atomtypes):
    """Write .def and .cif files for a randomly-generated porous material.

    Args:
        run_id (str): identification string for run.
        number_of_atomtypes (int): number of different chemical species used to
            populate the unit cell.

    Returns:
        material (sqlalchemy.orm.query.Query): database row for storing 
            simulation data specific to the material. See
            `htsohm/db/material.py` for more information.

        Atom-sites and unit-cell dimensions are stored in a .cif-file within
        RASPA's structures library: `$(raspa-dir)/structures/cif`. Force field
        parameters are stored within that material's directory in the RASPA
        library: `$(raspa-dir)/forcefield/(run_id)-(uuid)`. Partial charges,
        atomic radii, and more define each chemical species in 
        `pseudo_atoms.def`. Lennard-Jones type interactions (sigma, 
        epsilon-values) are defined in `force_field_mixing_rules.def`. These
        interactions can be overwritten in `force_field.def`, but by default
        no interactions are overwritten in this file.
 
    """
    lattice_limits          = config["lattice_constant_limits"]
    number_density_limits   = config["number_density_limits"]
    epsilon_limits          = config["epsilon_limits"]
    sigma_limits            = config["sigma_limits"]
    max_charge              = config["charge_limit"]
    elem_charge             = config["elemental_charge"]
    raspa2_dir              = config["raspa2_dir"]
    ff_dir   = os.path.join(raspa2_dir, 'share', 'raspa', 'forcefield')         # forcefield-files
    mat_dir  = os.path.join(raspa2_dir, 'share', 'raspa', 'structures', 'cif')  # .cif-files

    ########################################################################
    material = Material(run_id)
    material_name = run_id + '-' + material.uuid
    material.generation = 0

    def_dir = os.path.join(ff_dir, material_name)       # directory for material's force field
    os.mkdir(def_dir)
    force_field_file = os.path.join(def_dir, 'force_field.def')      # for overwriting LJ-params
    write_force_field(force_field_file)

    ########################################################################
    # define pseudo atom types by randomly-generating sigma and epsilon values
    #
    # NOTE :  assignment of partial charges is not supported in this version,
    #         and not necessary when modelling gases without dipole moments
    #         (such as methane). Charge assignment and mutation will appear
    #         in the next stable release.
    #
    atom_types = []
    for chemical_id in range(number_of_atomtypes):
        atom_types.append({
            "chemical-id" : "A_%s" % chemical_id,
            "charge"      : 0.,    # See NOTE above.
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

    return material

def closest_distance(x, y):
    """Finds closest distance between two points across periodic boundaries.

    Args:
        x (float): value in range [0,1].
        y (float): value in range [0,1].

    Returns:
        D (float): closest distance measured between `x` and `y`, with
            periodic boundaries at 0 and 1. For example, the disance between
            x = 0.2 and y = 0.8 would be 0.4, and not 0.6.

    """
    a = 1 - y + x
    b = abs(y - x)
    c = 1 - x + y
    return min(a, b, c)

def random_position(x_o, x_r, strength):
    """Produces a point along the closest path between two points.

    Args:
        x_o (float): value in range [0,1].
        x_r (float): value in range [0,1].
        strength (float): refers to mutation strength. Value determines
            fractional distance to `x_r` from `x_o`.

    Returns:
        xfrac (float): a value between `x_o` and `x_r` across the closest
            distance accross periodic boundaries at 0 and 1.

    """
    dx = closest_distance(x_o, x_r)
    if (x_o > x_r
            and (x_o - x_r) > 0.5):
        xfrac = round((x_o + strength * dx) % 1., 4)
    if (x_o < x_r
            and (x_r - x_o) > 0.5):
        xfrac = round((x_o - strength * dx) % 1., 4)
    if (x_o >= x_r
            and (x_o - x_r) < 0.5):
        xfrac = round(x_o - strength * dx, 4)
    if (x_o < x_r
            and (x_r - x_o) < 0.5):
        xfrac = round(x_o + strength * dx, 4)
    return xfrac

def write_child_definition_files(run_id, parent_id, generation, mutation_strength):
    """Modifies a "parent" material's definition files by perturbing each
    parameter by some factor, dictated by the `mutation_strength`.
    
    Args:
        run_id (str): identification string for run.
        parent_id (str): uuid, identifying parent material in database.
        generation (int): iteration count for overall bin-mutate-simulate routine.
        mutation_strength (float): perturbation factor [0, 1].

    Returns:
        new_material (sqlalchemy.orm.query.Query): database row for storing 
            simulation data specific to the material. See
            `htsohm/db/material.py` for more information.

    Todo:
        * Add methods for assigning and mutating charges.

    """
    parent = session.query(Material).get(str(parent_id))

    ########################################################################
    # load boundaries from config-file
    lattice_limits          = config["lattice_constant_limits"]
    number_density_limits   = config["number_density_limits"]
    epsilon_limits          = config["epsilon_limits"]
    sigma_limits            = config["sigma_limits"]
    raspa2_dir              = config["raspa2_dir"]
    ff_dir   = os.path.join(raspa2_dir, 'share', 'raspa', 'forcefield')         # forcefield-files
    mat_dir  = os.path.join(raspa2_dir, 'share', 'raspa', 'structures', 'cif')  # .cif-files

    ########################################################################
    # add row to database
    new_material = Material(run_id)
    new_material.parent_id = parent_id
    new_material.generation = generation

    ########################################################################
    # write force_field.def
    child_forcefield_directory = os.path.join(ff_dir, run_id + '-' + str(new_material.uuid))
    os.mkdir(child_forcefield_directory)
    force_field_file = os.path.join(child_forcefield_directory, 'force_field.def')
    write_force_field(force_field_file)                        # WRITE FORCE_FIELD.DEF
    print('  - force_field.def\tWRITTEN.')

    ########################################################################
    # copy pseudo_atoms.def
    parent_forcefield_directory = os.path.join(ff_dir, run_id + '-' + str(parent.uuid))
    parent_pseudo_file = os.path.join(parent_forcefield_directory, 'pseudo_atoms.def')
    shutil.copy(parent_pseudo_file, child_forcefield_directory)    # COPY PSEUDO_ATOMS.DEF
    print('  - pseudo_atoms.def\tWRITTEN.')

    ########################################################################
    # perturb LJ-parameters, write force_field_mixing_rules.def
    p_mix = os.path.join(parent_forcefield_directory, 'force_field_mixing_rules.def')
    n1, n2, old_epsilons, old_sigmas = np.genfromtxt(
        p_mix, unpack=True, skip_header=7, skip_footer=9
    )
    chemical_ids = np.genfromtxt(
        p_mix, unpack=True, skip_header=7, skip_footer=9, usecols=0, dtype=str
    )
    atom_types = []
    for i in range(len(chemical_ids)):
        epsilon = round( old_epsilons[i] + mutation_strength * (uniform(*epsilon_limits) -
        old_epsilons[i]), 4)
        sigma = round(
            old_sigmas[i] + mutation_strength * (uniform(*sigma_limits) - old_sigmas[i]
        ), 4)
        atom_type = {
            "chemical-id" : chemical_ids[i],
            "charge"      : 0.,    # See Docstring Todo.
            "epsilon"     : epsilon,
            "sigma"       : sigma}
        atom_types.append(atom_type)

        mix_file = os.path.join(child_forcefield_directory, 'force_field_mixing_rules.def')
        write_mixing_rules(mix_file, atom_types)
        print('  - force_field_mixing_rules.def\tWRITTEN.')

    ########################################################################
    # load values from parent cif-file
    p_cif = os.path.join(mat_dir, run_id + '-' + str(parent.uuid) + '.cif')
    n1, n2, old_x, old_y, old_z = np.genfromtxt(p_cif, unpack=True, skip_header=16)
    old_atom_types = np.genfromtxt(p_cif, usecols=0, dtype=str, skip_header=16)
    old_abc = np.genfromtxt(
        p_cif, unpack=True, skip_header=4, skip_footer=len(old_x) + 9, usecols = 1
    )
    ########################################################################
    # calculate new lattice constants
    old_a = old_abc[0]
    old_b = old_abc[1]
    old_c = old_abc[2]
    random_a = round(uniform(*lattice_limits), 4)
    random_b = round(uniform(*lattice_limits), 4)
    random_c = round(uniform(*lattice_limits), 4)
    a = round(old_a + mutation_strength * (random_a - old_a), 4)
    b = round(old_b + mutation_strength * (random_b - old_b), 4)
    c = round(old_c + mutation_strength * (random_c - old_c), 4)
    lattice_constants = {"a" : a, "b" : b, "c" : c}
    ########################################################################
    # calulate new number density, number of atoms
    old_ND = len(old_x) / (old_a * old_b * old_c)
    random_ND = round(uniform(*number_density_limits), 4)
    number_density = old_ND + mutation_strength * (random_ND - old_ND)
    number_of_atoms = int(number_density * a * b * c)
    ########################################################################
    # remove excess atom-sites, if any
    if number_of_atoms < len(old_x):
        difference = len(old_x) - number_of_atoms
        old_x = old_x[:-difference]
        old_y = old_y[:-difference]
        old_z = old_z[:-difference]
    ########################################################################
    # perturb atom-site positions
    atom_sites = []
    for i in range(len(old_x)):
        chemical_id = old_atom_types[i]
        x = random_position(old_x[i], random(), mutation_strength)
        y = random_position(old_y[i], random(), mutation_strength)
        z = random_position(old_z[i], random(), mutation_strength)
        atom_site = {
            "chemical-id" : chemical_id,
            "x-frac"      : x,
            "y-frac"      : y,
            "z-frac"      : z}
        atom_sites.append(atom_site)
    ########################################################################
    # add atom-sites, if needed
    if number_of_atoms > len(atom_sites):
        for new_sites in range(number_of_atoms - len(atom_sites)):
            atom_site = {
                "chemical-id" : choice(chemical_ids),
                "x-frac"      : round(random(), 4),
                "y-frac"      : round(random(), 4),
                "z-frac"      : round(random(), 4)}
            atom_sites.append(atom_site)

    cif_file = os.path.join(mat_dir, run_id + '-' + str(new_material.uuid) + '.cif')
    write_cif_file(cif_file, lattice_constants, atom_sites)

    return new_material

def write_cif_file(cif_file, lattice_constants, atom_sites):
    """Writes .cif file for structural information.

    Args:
        cif_file (str): path to .cif-file, for example:
            `$(raspa-dir)/structures/cif/(run_id)-(uuid).cif`.
        lattice_constants (dict): crystal lattice constants, for example:
            {"a" : (float),
             "b" : (float),
             "c" : (float)}
        atom_sites (list): dictionaries for each atom-site describing position
            and chemical species, for example:
            {"chemical-id" : chemical_id,
             "x-frac"      : x,
             "y-frac"      : y,
             "z-frac"      : z}

        Returns:
            Writes .cif file in RASPA's library,
            `$(raspa-dir)/structures/cif/(run_id)-(uuid).cif`.

    """
    with open(cif_file, "w") as file:
        file.write(
            "\nloop_\n" +
            "_symmetry_equiv_pos_as_xyz\n" +
            "  x,y,z\n" +
            "_cell_length_a\t%s\n" % lattice_constants["a"] +
            "_cell_length_b\t%s\n" % lattice_constants["b"] +
            "_cell_length_c\t%s\n" % lattice_constants["c"] +
            "_cell_angle_alpha\t90.0000\n" +
            "_cell_angle_beta\t90.0000\n" +
            "_cell_angle_gamma\t90.0000\n" +
            "loop_\n" +
            "_atom_site_label\n" +
            "_atom_site_type_symbol\n" +
            "_atom_site_fract_x\n" +
            "_atom_site_fract_y\n" +
            "_atom_site_fract_z\n")
        for atom_site in atom_sites:
            chemical  = atom_site["chemical-id"]
            x         = atom_site["x-frac"]
            y         = atom_site["y-frac"]
            z         = atom_site["z-frac"]
            file.write("%s\tC\t%s\t%s\t%s\n" % (chemical, x, y, z))

def write_mixing_rules(mix_file, atom_types):
    """Writes .def file for forcefield information.

    Args:
        mix_file (str): path to mixing-rules file, for example:
            `$(raspa-dir)/forcefield/(run_id)-(uuid)/force_field_mixing_rules.def`
        atom_types (list): dictionaries for each atom-type describing LJ-type
            interactions, for example:
            {"chemical-id" : chemical_ids[i],
             "charge"      : 0.,
             "epsilon"     : epsilon,
             "sigma"       : sigma}

    Returns:
        Writes file within RASPA's library,
        `$(raspa-dir)/forcefield/(run_id)-(uuid)/force_field_mixing_rules.def`

    """
    with open(mix_file, "w") as file:
        file.write(
            "# general rule for shifted vs truncated\n" +
            "shifted\n" +
            "# general rule tailcorrections\n" +
            "no\n" +
            "# number of defined interactions\n" +
            "%s\n" % (len(atom_types) + 8) +
            "# type interaction, parameters.    " +
            "IMPORTANT: define shortest matches first, so" +
            " that more specific ones overwrites these\n")
        for atom_type in atom_types:
            atom_id    = atom_type["chemical-id"]
            epsilon    = atom_type["epsilon"]
            sigma      = atom_type["sigma"]
            file.write(
                "%s\tlennard-jones\t%s\t%s\n" % (atom_id, epsilon, sigma))
        file.write(
            "N_n2\tlennard-jones\t36.0\t3.31\n" +
            "N_com\tnone\n" +
            "C_co2\tlennard-jones\t27.0\t2.80\n" +
            "O_co2\tlennard-jones\t79.0\t3.05\n" +
            "CH4_sp3\tlennard-jones\t158.5\t3.72\n" +
            "He\tlennard-jones\t10.9\t2.64\n" +
            "H_h2\tnone\n" +
            "H_com\tlennard-jones\t36.7\t2.958\n" +
            "# general mixing rule for Lennard-Jones\n" +

            "Lorentz-Berthelot")
def write_pseudo_atoms(psu_file, atom_types):
    """Writes .def file for chemical information.

    Args:
        psu_file (str): path to pseudo atoms definitions file, for example:
            `$(raspa-dir)/forcefield/(run_id)-(uuid)/pseudo_atoms.def`
        atom_types (list): dictionaries for each chemical species, including
            an identification string and charge, for example:
            interactions, for example:
            {"chemical-id" : chemical_ids[i],
             "charge"      : 0.,
             "epsilon"     : epsilon,
             "sigma"       : sigma}

    Returns:
        Writes file within RASPA's library,
        `$(raspa-dir)/forcefield/(run_id)-(uuid)/pseudo_atoms.def`

        NOTE: ALL CHARGES ARE 0. IN THIS VERSION.

    """
    with open(psu_file, "w") as file:
        file.write(
            "#number of pseudo atoms\n" +
            "%s\n" % (len(atom_types) + 8) +
            "#type\tprint\tas\tchem\toxidation\tmass\tcharge\tpolarization\tB-factor\tradii\t" +
                 "connectivity\tanisotropic\tanisotrop-type\ttinker-type\n")
        for atom_type in atom_types:
            atom_id   = atom_type["chemical-id"]
            charge    = atom_type["charge"]
            file.write(
                "A_%s\tyes\tC\tC\t0\t12.\t%s\t0\t0\t1\t1\t0\t0\tabsolute\t0\n" % (atom_id, charge))
        file.write(
            "N_n2\tyes\tN\tN\t0\t14.00674\t-0.4048\t0.0\t1.0\t0.7\t0\t0\trelative\t0\n" +
            "N_com\tno\tN\t-\t0\t0.0\t0.8096\t0.0\t1.0\t0.7\t0\t0\trelative\t0\n" +
            "C_co2\tyes\tC\tC\t0\t12.0\t0.70\t0.0\t1.0\t0.720\t0\t0\trelative\t0\n" +
            "O_co2\tyes\tO\tO\t0\t15.9994\t-0.35\t0.0\t1.0\t0.68\t0\t0\trelative\t0\n" +
            "CH4_sp3\tyes\tC\tC\t0\t16.04246\t0.0\t0.0\t1.0\t1.00\t0\t0\trelative\t0\n" +
            "He\tyes\tHe\tHe\t0\t4.002602\t0.0\t0.0\t1.0\t1.0\t0\t0\trelative\t0\n" +
            "H_h2\tyes\tH\tH\t0\t1.00794\t0.468\t0.0\t1.0\t0.7\t0\t0\trelative\t0\n" +
            "H_com\tno\tH\tH\t0\t0.0\t-0.936\t0.0\t1.0\t0.7\t0\t0\trelative\t0\n")

def write_force_field(for_file):
    """Writes .def file to overwrite LJ-type interactions.

    Args:
        for_file (str): path to .def-file, for example:
            `$(raspa-dir)/forcefield/(run_id)-(uuid)/force_field.def`

    Returns:
        Writes file within RASPA's library.
        `$(raspa-dir)/forcefield/(run_id)-(uuid)/force_field.def`

        NOTE: NO INTERACTIONS ARE OVERWRITTEN BY DEFAULT.

    """
    with open(for_file, "w") as file:
        file.write(
            "# rules to overwrite\n" +
            "0\n" +
            "# number of defined interactions\n" +
            "0\n" +
            "# mixing rules to overwrite\n" +
            "0")
