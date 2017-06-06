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
import htsohm
from htsohm import config
from htsohm.pseudo_material import PseudoMaterial
from htsohm.db import session, Material, Structure, LennardJones, AtomSites

def random_number_density(number_density_limits, material_structure): #lattice_constants):
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
    v = material_structure.volume
    min_atoms = int(min_ND * v)
    max_atoms = int(max_ND * v)
    if min_atoms < 2:
        min_atoms = int(2)
    atoms = randrange(min_atoms, max_atoms + 1, 1)
    return atoms

def generate_material(run_id, number_of_atomtypes):
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
    lattice_limits          = config["lattice_constant_limits"]
    number_density_limits   = config["number_density_limits"]
    epsilon_limits          = config["epsilon_limits"]
    sigma_limits            = config["sigma_limits"]
    max_charge              = config["charge_limit"]
    elem_charge             = config["elemental_charge"]

    ########################################################################
    material = Material(run_id)
    material.generation = 0

    structure = material.structure
    structure.lattice_constant_a = uniform(*lattice_limits)
    structure.lattice_constant_b = uniform(*lattice_limits)
    structure.lattice_constant_c = uniform(*lattice_limits)

    for chemical_id in range(number_of_atomtypes):
        structure.lennard_jones.append(LennardJones(
                    chemical_id = 'A_{}'.format(chemical_id),
                    sigma = uniform(*sigma_limits),
                    epsilon = uniform(*epsilon_limits)))

    atom_sites = structure.atom_sites
    for i in range(random_number_density(number_density_limits, structure)):
        atom_sites.append(AtomSites(
            chemical_id = 'A_{}'.format(choice(range(number_of_atomtypes))),
            x_frac = random(), y_frac = random(), z_frac = random()))

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
        xfrac = (x_o + strength * dx) % 1.
    if (x_o < x_r
            and (x_r - x_o) > 0.5):
        xfrac = (x_o - strength * dx) % 1.
    if (x_o >= x_r
            and (x_o - x_r) < 0.5):
        xfrac = x_o - strength * dx
    if (x_o < x_r
            and (x_r - x_o) < 0.5):
        xfrac = x_o + strength * dx
    return xfrac

def mutate_material(parent_material, mutation_strength, generation):
    """Modifies a "parent" pseudomaterial structure by perturbing each
    parameter by some factor, dictated by the `mutation_strength`.
    
    Args:
        parent_material : record for parent pseudomaterial in 'materials' table.
        mutation_strength (float): perturbation factor [0, 1].
        generation (int): iteration count for overall bin-mutate-simulate routine.

    Returns:

    Todo:
        * Add methods for assigning and mutating charges.

    """
    
    ########################################################################
    # load boundaries from config-file
    lattice_limits          = config["lattice_constant_limits"]
    number_density_limits   = config["number_density_limits"]

    child_material = Material(parent_material.run_id)
    child_material.generation = generation
    child_material.parent_id = parent_material.id
 
    # perturb lennard-jones parameters
    for atom_type in parent_material.structure.lennard_jones:
        child_material.structure.lennard_jones.append(LennardJones(
            chemical_id = atom_type.chemical_id,
            sigma = atom_type.sigma + mutation_strength * (uniform(*config['sigma_limits']) - atom_type.sigma),
            epsilon = atom_type.epsilon + mutation_strength * (uniform(*config['epsilon_limits']) - atom_type.epsilon)))

    # perturb lattice constants
    child_material.structure.lattice_constant_a = parent_material.structure.lattice_constant_a \
            + mutation_strength * (uniform(*lattice_limits) - parent_material.structure.lattice_constant_a)
    child_material.structure.lattice_constant_b = parent_material.structure.lattice_constant_b \
            + mutation_strength * (uniform(*lattice_limits) - parent_material.structure.lattice_constant_b)
    child_material.structure.lattice_constant_c = parent_material.structure.lattice_constant_c \
            + mutation_strength * (uniform(*lattice_limits) - parent_material.structure.lattice_constant_c)

    # perturb number density/number of atom-sites
    parent_ND = len(parent_material.structure.atom_sites) / parent_material.structure.volume
    child_ND = parent_ND + mutation_strength * (uniform(*number_density_limits) - parent_ND)
    number_of_atoms = int(child_ND * child_material.structure.volume)

    # remove atom-sites, if necessary
    child_material.structure.atom_sites = np.random.choice(
        parent_material.structure.atom_sites,
        min(number_of_atoms, len(parent_material.structure.atom_sites)),
        replace = False).tolist()

    # perturb atom-site positions
    for atom_site in child_material.structure.atom_sites:
        atom_site.x_frac = random_position(atom_site.x_frac, random(), mutation_strength)
        atom_site.y_frac = random_position(atom_site.y_frac, random(), mutation_strength)
        atom_site.z_frac = random_position(atom_site.z_frac, random(), mutation_strength)

    # add atom-sites, if needed
    if number_of_atoms > len(parent_material.structure.atom_sites):
        for new_sites in range(number_of_atoms - len(parent_material.structure.atom_sites)):
            parent_material.structure.atom_sites.append(AtomSites(
                chemical_id = 'A_{}'.format(choice(
                    range(len(parent_material.structure.lennard_jones)))),
                x_frac = random(), y_frac = random(), z_frac = random()))

    return child_material

def write_cif_file(material, simulation_path):
    """Writes .cif file for structural information.

    Args:


    """
    file_name = os.path.join(simulation_path, '%s.cif' % material.uuid)
    with open(file_name, "w") as cif_file:
        cif_file.write(
            "\nloop_\n" +
            "_symmetry_equiv_pos_as_xyz\n" +
            "  x,y,z\n" +
            "_cell_length_a\t{}\n".format(round(material.structure.lattice_constant_a)) +
            "_cell_length_b\t{}\n".format(round(material.structure.lattice_constant_b)) +
            "_cell_length_c\t{}\n".format(round(material.structure.lattice_constant_c)) +
            "_cell_angle_alpha  90.0000\n" +
            "_cell_angle_beta   90.0000\n" +
            "_cell_angle_gamma  90.0000\n" +
            "loop_\n" +
            "_atom_site_label\n" +
            "_atom_site_type_symbol\n" +
            "_atom_site_fract_x\n" +
            "_atom_site_fract_y\n" +
            "_atom_site_fract_z\n"
        )
        for atom_site in material.structure.atom_sites:
            cif_file.write(
            "{0:5} C {1:4f} {2:4f} {3:4f}\n".format(
                atom_site.chemical_id,
                round(atom_site.x_frac, 4),
                round(atom_site.y_frac, 4),
                round(atom_site.z_frac, 4)
            ))

def write_mixing_rules(material, simulation_path):
    """Writes .def file for forcefield information.

    """
    adsorbate_LJ_atoms = [
            ['N_n2',    36.0,       3.31],
            ['C_co2',   27.0,       2.80],
            ['O_co2',   79.0,       3.05],
            ['CH4_sp3', 158.5,      3.72],
            ['He',      10.9,       2.64],
            ['H_com',   36.7,       2.958],
            ['Kr',      167.06,     3.924],
            ['Xe',      110.704,    3.690]
    ]
 
    adsorbate_none_atoms = ['N_com', 'H_h2']

    file_name = os.path.join(simulation_path, 'force_field_mixing_rules.def')
    with open(file_name, "w") as mixing_rules_file:
        mixing_rules_file.write(
            "# general rule for shifted vs truncated\n" +
            "shifted\n" +
            "# general rule tailcorrections\n" +
            "no\n" +
            "# number of defined interactions\n" +
            "{}\n".format(len(material.structure.lennard_jones) + 10) +
            "# type interaction, parameters.    " +
            "IMPORTANT: define shortest matches first, so" +
            " that more specific ones overwrites these\n"
        )
        for lennard_jones in material.structure.lennard_jones:
            mixing_rules_file.write(
                "{0:12} lennard-jones {1:8f} {2:8f}\n".format(
                    lennard_jones.chemical_id,
                    round(lennard_jones.epsilon, 4),
                    round(lennard_jones.sigma, 4)
                )
            )
        for at in adsorbate_LJ_atoms:
            mixing_rules_file.write(
                "{0:12} lennard-jones {1:8f} {2:8f}\n".format(at[0], at[1], at[2])
            )
        for at in adsorbate_none_atoms:
            mixing_rules_file.write(
                "{0:12} none\n".format(at)
            )
        mixing_rules_file.write(
            "# general mixing rule for Lennard-Jones\n" +
            "Lorentz-Berthelot")

def write_pseudo_atoms(material, simulation_path):
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
    temporary_charge = 0.

    file_name = os.path.join(simulation_path, 'pseudo_atoms.def')
    with open(file_name, "w") as pseudo_atoms_file:
        pseudo_atoms_file.write(
            "#number of pseudo atoms\n" +
            "%s\n" % (len(material.structure.lennard_jones) + 10) +
            "#type  print   as  chem    oxidation   mass    charge  polarization    B-factor    radii   " +
                 "connectivity  anisotropic anisotrop-type  tinker-type\n")
        for atom_type in material.structure.lennard_jones:
            pseudo_atoms_file.write(
                "{0:7}  yes  C   C   0   12.0       {0:8}  0.0  1.0  1.0    0  0  absolute  0\n".format(
                    atom_type.chemical_id, 0.0
                )
            )
        pseudo_atoms_file.write(
            "N_n2     yes  N   N   0   14.00674   -0.4048   0.0  1.0  0.7    0  0  relative  0\n" +
            "N_com    no   N   -   0    0.0        0.8096   0.0  1.0  0.7    0  0  relative  0\n" +
            "C_co2    yes  C   C   0   12.0        0.70     0.0  1.0  0.720  0  0  relative  0\n" +
            "O_co2    yes  O   O   0   15.9994    -0.35     0.0  1.0  0.68   0  0  relative  0\n" +
            "CH4_sp3  yes  C   C   0   16.04246    0.0      0.0  1.0  1.00   0  0  relative  0\n" +
            "He       yes  He  He  0    4.002602   0.0      0.0  1.0  1.0    0  0  relative  0\n" +
            "H_h2     yes  H   H   0    1.00794    0.468    0.0  1.0  0.7    0  0  relative  0\n" +
            "H_com    no   H   H   0    0.0        0.936    0.0  1.0  0.7    0  0  relative  0\n" +
            "Xe       yes  Xe  Xe  0  131.293      0.0      0.0  1.0  2.459  0  0  relative  0\n" +
            "Kr       yes  Kr  Kr  0   83.798      0.0      0.0  1.0  2.27   0  0  relative  0\n"
        )

def write_force_field(simulation_path):
    """Writes .def file to overwrite LJ-type interactions.

    Args:
        file_name (str): path to .def-file, for example:
            `$(raspa-dir)/forcefield/(run_id)-(uuid)/force_field.def`

    Writes file within RASPA's library:
        `$(raspa-dir)/forcefield/(run_id)-(uuid)/force_field.def`

    NOTE: NO INTERACTIONS ARE OVERWRITTEN BY DEFAULT.

    """
    file_name = os.path.join(simulation_path, 'force_field.def')
    with open(file_name, "w") as force_field_file:
        force_field_file.write(
            "# rules to overwrite\n" +
            "0\n" +
            "# number of defined interactions\n" +
            "0\n" +
            "# mixing rules to overwrite\n" +
            "0")
