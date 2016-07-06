# stanard imports
import os
from math import sqrt

# related third party imports
from datetime import datetime
import yaml
from sqlalchemy import func

# local application/library specific imports
from htsohm.runDB_declarative import Base, Material, session

def write_config_file(children_per_generation, number_of_atomtypes, strength_0,
    number_of_bins, max_generations, dummy_test_trials, acceptance_value):
    """Writes run-specific parameters to /config/<run_id>.yaml.

    This function writes run-specific parameters to a configuration file that is loaded/written
    at different stages throughout the overall HTSOHM routine.
    """
    run_id = datetime.now().isoformat()
    wd = os.environ['HTSOHM_DIR']      # specify working directory
    config_file = os.path.join(wd, 'config', run_id + '_conf.yaml')

    run_config = {
        "run-id" : run_id,
        "children-per-generation" : children_per_generation,
        "number-of-atom-types" : number_of_atomtypes,
        "initial-mutation-strength" : strength_0,
        "number-of-bins" : number_of_bins,
        "max-number-of-generations" : max_generations,
        "dummy-test-trials" : dummy_test_trials,
        "acceptance-value" : acceptance_value
    }

    with open(config_file, "w") as file:
        yaml.dump(run_config, file, default_flow_style=False)

    return run_config

def update_config_file(run_id):
    """ Write material-parameters to run-configuration file.

    The parameters written by this function define the limits for different values written to the
    structure and forcefield definition files for RASPA. Among the limits defined here are crystal
    lattice constants, number density, partial atomic charges, and Lennard-Jones parameters (sigma
    and epsilon).
    """
    wd = os.environ['HTSOHM_DIR']      # specify $HTSOHM_DIR as working directory
    config_file = os.path.join(wd, 'config', run_id + '_conf.yaml')
    with open(config_file) as file:
        run_config = yaml.load(file)

    run_config.update({
        "number-density-limits"     : [0.0000149, 0.02122],
        "lattice-constant-limits"   : [25.6, 51.2],
        "epsilon-limits"            : [1.258, 513.264],
        "sigma-limits"              : [1.052, 6.549],
        "charge-limit"              : 0.,
        "elemental-charge"          : 0.0001
    })

    with open(config_file, "w") as file:
        yaml.dump(run_config, file, default_flow_style=False)

    return run_config

def read_config_file(run_id):
    wd = os.environ['HTSOHM_DIR']
    config_file = os.path.join(wd, 'config', run_id + '_conf.yaml')
    with open(config_file) as file:
        config = yaml.load(file)
    return config

def evaluate_convergence(run_id):
    '''Counts number of materials in each bin and returns variance of these counts.'''
    bin_counts = session \
        .query(func.count(Material.id)) \
        .filter(Material.run_id == run_id, Material.write_check == 'done') \
        .group_by(
            Material.methane_loading_bin, Material.surface_area_bin, Material.void_fraction_bin
        ).all()
    bin_counts = [i[0] for i in bin_counts]    # convert SQLAlchemy result to list
    variance = sqrt( sum([(i - (sum(bin_counts) / len(bin_counts)))**2 for i in bin_counts]) / len(bin_counts))
    return variance

def save_convergence(run_id, generation, variance):
    wd = os.environ['HTSOHM_DIR']      # specify working directory
    log_file = os.path.join(wd, 'config', run_id + '_log.yaml')
    data = {
        "generation_%s" % generation  : variance
    }
    if generation == 0:
        with open(log_file, "w") as file:
            yaml.dump(data, file, default_flow_style=False)
    else:
        with open(log_file, "a") as file:
            yaml.dump(data, file, default_flow_style=False)

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
