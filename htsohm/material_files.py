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
from htsohm import config, PseudoMaterial
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
    atoms = randrange(min_atoms, max_atoms + 1, 1)
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

    Atom-sites and unit-cell dimensions are stored in a .cif-file within RASPA's
    structures library: `$(raspa-dir)/structures/cif`. Force field parameters
    are stored within that material's directory in the RASPA library: 
    `$(raspa-dir)/forcefield/(run_id)-(uuid)`. Partial charges, atomic radii,
    and more define each chemical species in `pseudo_atoms.def`. Lennard-Jones
    type interactions (sigma, epsilon-values) are defined in 
    `force_field_mixing_rules.def`. These interactions can be overwritten in 
    `force_field.def`, but by default no interactions are overwritten in this file.
 
    """
    lattice_limits          = config["lattice_constant_limits"]
    number_density_limits   = config["number_density_limits"]
    epsilon_limits          = config["epsilon_limits"]
    sigma_limits            = config["sigma_limits"]
    max_charge              = config["charge_limit"]
    elem_charge             = config["elemental_charge"]

    htsohm_dir = os.path.dirname(os.path.dirname(htsohm.__file__))
    run_dir = os.path.join(htsohm_dir, run_id)
    material_dir = os.path.join(run_dir, 'pseudo_materials')
    if not os.path.exists(material_dir):
        os.mkdir(material_dir)

    ########################################################################
    db_row = Material(run_id)
    db_row.generation = 0
    material = PsuedoMaterial(db_row.uuid)

    material.atom_types = []
    for chemical_id in range(number_of_atomtypes):
        material.atom_types.append({
            "chemical-id" : "A_%s" % chemical_id,
            "charge"      : 0.,    # See NOTE above.
            "epsilon"     : round(uniform(*epsilon_limits), 4),
            "sigma"       : round(uniform(*sigma_limits), 4)
        })

    material.lattice_constants = {}
    for i in ['a', 'b', 'c']:
        material.lattice_constants[i] = round(uniform(*lattice_limits), 4)

    material.number_of_atoms   = random_number_density(
        number_density_limits, lattice_constants)

    material.atom_sites = []
    for atom in range(number_of_atoms):
        atom_site = {"chemical-id" : choice(atom_types)["chemical-id"]}
        for i in ['x-frac', 'y-frac', 'z-frac']:
            atom_site[i] = round(random(), 4)
        material.atom_sites.append(atom_site)

    material_file = os.path.join(material_dir, '%s.yaml' % uuid)
    with open(material_file) as dump_file:
        yaml.dump(material, dump_file) 

    return db_row

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
    
    ########################################################################
    # load boundaries from config-file
    lattice_limits          = config["lattice_constant_limits"]
    number_density_limits   = config["number_density_limits"]
    epsilon_limits          = config["epsilon_limits"]
    sigma_limits            = config["sigma_limits"]

    htsohm_dir = os.path.dirname(os.path.dirname(htsohm.__file__))
    run_dir = os.path.join(htsohm_dir, run_id)
    material_dir = os.path.join(run_dir, 'pseudo_materials')

    parent_row = session.query(Material).get(str(parent_id))
    parent_yaml = os.path.join(material_dir, '%s.yaml' % parent_row.uuid)
    with open(parent_yaml) as load_file:
        parent_material = yaml.load(load_file)

    ########################################################################
    # add row to database
    child_row = Material(run_id)
    child_row.parent_id = parent_id
    child_row.generation = generation
    child_material = PseudoMaterial(new_db_row.uuid)

    ########################################################################
    # perturb LJ-parameters
    child_material.atom_types = []
    for atom_type in parent_material.atom_types:    
        new_atom_type = {'chemical-id' : atom_type['chemical-id']}
        new_atom_type['epsilon'] = (round(atom_type['epsilon'] +
                mutation_strength * (uniform(*epsilon_limits) -
                atom_type['epsilon']), 4))
        new_atom_type['sigma'] = (round(atom_type['sigma'] +
                mutation_strength * (uniform(*epsilon_limits) -
                atom_type['epsilon']), 4))
        child_material.atom_types.append(new_atom_type)

    ########################################################################
    # calculate new lattice constants
    child_material.lattice_constants = {}
    for i in ['a', 'b', 'c']:
        old_x = parent_material.lattice_constants[i]
        random_x = uniform(*lattice_limits)
        new_x = round(old_x + mutation_strength * (random_x - old_x), 4)
        child_material.lattice_constants[i] = new_x

    ########################################################################
    # calulate new number density, number of atoms
    old_LCs = parent_material.lattice_constants
    old_vol = old_LCs['a'] * old_LCs['b'] * old_LCs['c']
    old_number_density = len(parent_material.atom_sites) / old_vol
    random_number_density = uniform(*number_density_limits)
    new_number_density = (old_number_density + mutation_strength * (
            random_number_density - old_number_density))
    new_LCs = child_material.lattice_constants
    new_vol = new_LCs['a'] * new_LCs['b'] * new_LCs['c']
    number_of_atom_sites = int(new_number_density * new_vol)

    ########################################################################
    # remove excess atom-sites, if any
    if number_of_atom_sites < len(parent_material.atom_sites):
        parent_material.atom_sites = parent_material.atom_sites[
                :number_of_atom_sites]
    ########################################################################
    # perturb atom-site positions
    child_material.atom_sites = []
    for atom_site in parent_material.atom_sites:
        new_atom_site = {'chemical-id' : atom_site}
        for i in ['x-frac', 'y-frac', 'z-frac']
            new_atom_site[i] = random_position(
                    atom_site[i], random(), mutation_strength)
        child_material.atom_sites.append(new_atom_site)

    ########################################################################
    # add atom-sites, if needed
    if number_of_atom_sites > len(child_material.atom_sites):
        for new_sites in range(number_of_atoms - len(atom_sites)):
            new_atom_site = {
                    'chemical-id' : choice(
                        child_material.atom_types)['chemical-id']}
            for i in ['x-frac', 'y-frac', 'z-frac']:
                new_atom_site[i] = round(random(), 4)
            atom_sites.append(new_atom_site)

    material_file = os.path.join(material_dir, '%s.yaml' % uuid)
    with open(material_file) as dump_file:
        yaml.dump(child_material, dump_file) 

    return child_row

def write_cif_file(file_name, lattice_constants, atom_sites):
    """Writes .cif file for structural information.

    Args:
        file_name (str): path to .cif-file, for example:
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

    Writes .cif file in RASPA's library:
        `$(raspa-dir)/structures/cif/(run_id)-(uuid).cif`

    """
    with open(file_name, "w") as cif_file:
        cif_file.write(
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
            cif_file.write("%s\tC\t%s\t%s\t%s\n" % (chemical, x, y, z))

def write_mixing_rules(file_name, atom_types):
    """Writes .def file for forcefield information.

    Args:
        file_name (str): path to mixing-rules file, for example:
            `$(raspa-dir)/forcefield/(run_id)-(uuid)/force_field_mixing_rules.def`
        atom_types (list): dictionaries for each atom-type describing LJ-type
            interactions, for example:
            {"chemical-id" : chemical_ids[i],
             "charge"      : 0.,
             "epsilon"     : epsilon,
             "sigma"       : sigma}

    Writes file within RASPA's library:
        `$(raspa-dir)/forcefield/(run_id)-(uuid)/force_field_mixing_rules.def`

    """
    with open(file_name, "w") as mixing_rules_file:
        mixing_rules__file.write(
            "# general rule for shifted vs truncated\n" +
            "shifted\n" +
            "# general rule tailcorrections\n" +
            "no\n" +
            "# number of defined interactions\n" +
            "%s\n" % (len(atom_types) + 10) +
            "# type interaction, parameters.    " +
            "IMPORTANT: define shortest matches first, so" +
            " that more specific ones overwrites these\n")
        for atom_type in atom_types:
            atom_id    = atom_type["chemical-id"]
            epsilon    = atom_type["epsilon"]
            sigma      = atom_type["sigma"]
            file.write(
                "%s\tlennard-jones\t%s\t%s\n" % (atom_id, epsilon, sigma))
        mixing_rules_file.write(
            "N_n2\tlennard-jones\t36.0\t3.31\n" +
            "N_com\tnone\n" +
            "C_co2\tlennard-jones\t27.0\t2.80\n" +
            "O_co2\tlennard-jones\t79.0\t3.05\n" +
            "CH4_sp3\tlennard-jones\t158.5\t3.72\n" +
            "He\tlennard-jones\t10.9\t2.64\n" +
            "H_h2\tnone\n" +
            "H_com\tlennard-jones\t36.7\t2.958\n" +
            "Kr\tlennard-jones\t167.062\t3.924\n" +
            "Xe\tlennard-jones\t110.704\t3.690\n" +
            "# general mixing rule for Lennard-Jones\n" +
            "Lorentz-Berthelot")

def write_pseudo_atoms(file_name, atom_types):
    """Writes .def file for chemical information.

    Args:
        file_name (str): path to pseudo atoms definitions file, for example:
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
    with open(file_name, "w") as pseudo_atoms_file:
        pseudo_atoms_file.write(
            "#number of pseudo atoms\n" +
            "%s\n" % (len(atom_types) + 10) +
            "#type\tprint\tas\tchem\toxidation\tmass\tcharge\tpolarization\tB-factor\tradii\t" +
                 "connectivity\tanisotropic\tanisotrop-type\ttinker-type\n")
        for atom_type in atom_types:
            atom_id   = atom_type["chemical-id"]
            charge    = atom_type["charge"]
            pseudo_atoms_file.write(
                "A_%s\tyes\tC\tC\t0\t12.\t%s\t0\t0\t1\t1\t0\t0\tabsolute\t0\n" % (atom_id, charge))
        pseudo_atoms_file.write(
            "N_n2\tyes\tN\tN\t0\t14.00674\t-0.4048\t0.0\t1.0\t0.7\t0\t0\trelative\t0\n" +
            "N_com\tno\tN\t-\t0\t0.0\t0.8096\t0.0\t1.0\t0.7\t0\t0\trelative\t0\n" +
            "C_co2\tyes\tC\tC\t0\t12.0\t0.70\t0.0\t1.0\t0.720\t0\t0\trelative\t0\n" +
            "O_co2\tyes\tO\tO\t0\t15.9994\t-0.35\t0.0\t1.0\t0.68\t0\t0\trelative\t0\n" +
            "CH4_sp3\tyes\tC\tC\t0\t16.04246\t0.0\t0.0\t1.0\t1.00\t0\t0\trelative\t0\n" +
            "He\tyes\tHe\tHe\t0\t4.002602\t0.0\t0.0\t1.0\t1.0\t0\t0\trelative\t0\n" +
            "H_h2\tyes\tH\tH\t0\t1.00794\t0.468\t0.0\t1.0\t0.7\t0\t0\trelative\t0\n" +
            "H_com\tno\tH\tH\t0\t0.0\t-0.936\t0.0\t1.0\t0.7\t0\t0\trelative\t0\n" +
            "Xe\tyes\tXe\tXe\t0\t131.293\t0.0\t0.0\t1.0\t2.459\t0\t0\trelative\t0\n" +
            "Kr\tyes\tKr\tKr\t0\t83.798\t0.0\t0.0\t1.0\t2.27\t0\t0\trelative\t0\n"
        )

def write_force_field(file_name):
    """Writes .def file to overwrite LJ-type interactions.

    Args:
        file_name (str): path to .def-file, for example:
            `$(raspa-dir)/forcefield/(run_id)-(uuid)/force_field.def`

    Writes file within RASPA's library:
        `$(raspa-dir)/forcefield/(run_id)-(uuid)/force_field.def`

    NOTE: NO INTERACTIONS ARE OVERWRITTEN BY DEFAULT.

    """
    with open(file_name, "w") as force_field_file:
        force_field_file.write(
            "# rules to overwrite\n" +
            "0\n" +
            "# number of defined interactions\n" +
            "0\n" +
            "# mixing rules to overwrite\n" +
            "0")
