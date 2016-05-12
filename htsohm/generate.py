import os

from random import choice, random, randrange, randint
from functools import reduce
from math import fsum
import numpy as np
from math import floor
import yaml

def gcd(a, b):
    while b:
        a, b = b, a % b
    return a

def lcm(a, b):
    return a * b // gcd(a, b)

def lcmm(*args):
    return reduce(lcm, args)

def write_material_config(run_id):
    wd = os.environ['HTSOHM_DIR']      # specify $HTSOHM_DIR as working directory
    config_file = os.path.join(wd, 'config', run_id + '.yaml')
    material_config = {
        "number-density-limits"     : [0.000013907, 0.084086],
        "lattice-constant-limits"   : [13.098, 52.392],
        "epsilon-limits"            : [1.258, 513.264],
        "sigma-limits"              : [1.052, 6.549],
        "charge-limit"              : 0.,
        "elemental-charge"          : 0.0001}
    with open(config_file, "a") as file:
        yaml.dump(material_config, file, default_flow_style=False)
    return material_config

def random_value(limits):
    x_min = limits[0]
    x_max = limits[1]
    x = round(random() * (x_max - x_min) + x_min, 4)
    return x

def random_number_density(number_density_limits, lattice_constants):
    max_number_density = number_density_limits[1]
    a = lattice_constants["a"]
    b = lattice_constants["b"]
    c = lattice_constants["c"]
    max_number_of_atoms = int(max_number_density * a * b * c)
    number_of_atoms = randrange(2, max_number_of_atoms, 1)
    number_density = round(number_of_atoms / (a * b * c))
    return number_of_atoms, number_density

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
            file.write("%s\t%s\t%s\t%s\n" % (chemical, x, y, z))

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
                "A_%s\tlennard-jones\t%s\t%s\n" % (atom_id, epsilon, sigma))
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

def write_seed_definition_files(run_id, number_of_materials, number_of_atomtypes):
#    run_id                  = run_config["run-id"]
#    number_of_materials     = run_config["children-per-generation"]
#    number_of_atomtypes     = run_config["number-of-atom-types"]
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
    for material in range(number_of_materials):           # each iteration creates a new material
        material_id = run_id + '-' + str(material)

        def_dir = os.path.join(ff_dir, material_id)       # directory for material's force field
        os.mkdir(def_dir)
        mix_file = os.path.join(def_dir, 'force_field_mixing_rules.def') # LJ-parameters
        for_file = os.path.join(def_dir, 'force_field.def')              # to overwrite LJ
        psu_file = os.path.join(def_dir, 'pseudo_atoms.def')             # define atom-types
        write_force_field(for_file)

        cif_file = os.path.join(mat_dir, material_id + ".cif")              # structure file
        lattice_constants = {"a" : random_value(lattice_limits),
                             "b" : random_value(lattice_limits),
                             "c" : random_value(lattice_limits)}
        number_of_atoms, number_density = random_number_density(number_density_limits, 
            lattice_constants)

        atom_types = []
        chemical_id = 0
        for atom_type in range(number_of_atomtypes):
            epsilon = random_value(epsilon_limits)
            sigma = random_value(sigma_limits)
            charge = 0.
            atom_type = {
                "chemical-id" : chemical_id,
                "charge"      : charge,
                "epsilon"     : epsilon,
                "sigma"       : sigma}
            atom_types.append(atom_type)
            chemical_id += 1
        write_mixing_rules(mix_file, atom_types)
        write_pseudo_atoms(psu_file, atom_types)

        atom_sites = []
        for atom_site in range(number_of_atoms):
            atom_type = choice(atom_types)
            chemical_id = atom_type["chemical-id"]
            x_fraction = round(random(), 4)
            y_fraction = round(random(), 4)
            z_fraction = round(random(), 4)
            atom_site = {
                "chemical-id" : chemical_id,
                "x-frac"      : x_fraction,
                "y-frac"      : y_fraction,
                "z-frac"      : z_fraction}
            atom_sites.append(atom_site)
        write_cif_file(cif_file, lattice_constants, atom_sites)
        material = {"atom-types" : atom_types, "atom-sites" : atom_sites}
        materials.append(material)
    return materials
