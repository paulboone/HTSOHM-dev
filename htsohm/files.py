import os
import yaml

def default_configuration():
    return {
        'override_restart_errors': False,
        'keep_configs': False,
        'output_dir': os.getcwd(),
        'verbose': False,
        'void_fraction_subtype': 'raspa',
        'load_restart_path': False,
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
         config.update(yaml.load(config_file))

    enforce_config_ok(config)

    return config

def enforce_config_ok(config):
    assert config['void_fraction_subtype'] in ["raspa", "geo", "zeo"]
