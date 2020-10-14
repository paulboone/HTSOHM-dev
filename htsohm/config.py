import time
import os
import yaml

def default_configuration():
    return {
        'override_restart_errors': False,
        'keep_configs': False,
        'output_dir': os.getcwd(),
        'void_fraction_subtype': 'raspa',
        'num_processes': 1,
        'initial_points_random_seed': int(time.time()),
        'structure_parameters': {
            'minimum_site_distance': 0.0
        }
    }

def load_config_file(path):
    """Loads the config file.

    Args:
        path (str): config path.

    Returns:
        config (dict): parameters specified in config.
    """
    config = default_configuration()
    with open(path) as config_file:
         config.update(yaml.safe_load(config_file))

    enforce_config_ok(config)

    return config

def enforce_config_ok(config):
    assert config['void_fraction_subtype'] in ["raspa", "geo", "zeo"]
    assert config['selector_type'] in ["simplices-or-hull", "density-bin", "neighbor-bin",
                                        "best", "specific", "random"]
    assert config['structure_parameters']['lattice_cubic'] == True
    valid_perturbations = {"num_atoms", "atom_type_assignments", "atom_types", "lattice", "atom_sites", "charges"}
    assert set(config['structure_parameters']['perturb']) <= valid_perturbations
