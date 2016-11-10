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
    """Find the `closest` distance over a periodic boundary.
    """
    a = 1 - y + x
    b = abs(y - x)
    c = 1 - x + y
    return min(a, b, c)

def random_position(x_o, x_r, strength):
    """Given two values, return some random point between them.
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
    """Mutate a parent-material's definition files to create new, child-material.
    At this point, a generation of materials has been initialized with parent-material IDs (primary
    keys). This function loads the necessary parameters from the selected parent's definition
    files, perturbs all of these values, and then writes them to new definition-files, creating
    a new material. The degree to which each value is perturbed is controlled by the mutation-
    strength-parameter."""

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
    #
    # NOTE :  assignment of partial charges is not supported in this version,
    #         and not necessary when modelling gases without dipole moments
    #         (such as methane). Charge assignment and mutation will appear
    #         in the next stable release.
    #
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
            "charge"      : 0.,    # See NOTE above.
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
    with open(for_file, "w") as file:
        file.write(
            "# rules to overwrite\n" +
            "0\n" +
            "# number of defined interactions\n" +
            "0\n" +
            "# mixing rules to overwrite\n" +
            "0")
