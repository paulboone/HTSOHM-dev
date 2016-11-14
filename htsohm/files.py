# standard imports
import os

# related third party imports
import yaml

# local application/library specific imports
import htsohm

def load_config_file(file_name):
    """Reads input file.

    Args:
        file_name (str): auto-created config filename (ex. htsohm.sample.yaml).

    Returns:
        config (dict): parameters specified in config.

    """
    with open(file_name) as file:
         config = yaml.load(file)
    return config

def save_convergence_file(run_id, generation, variance):
    htsohm_dir = os.path.dirname(os.path.dirname(htsohm.__file__))
    log_file = os.path.join(htsohm_dir, 'config', run_id + '_log.yaml')
    data = {
        "generation_%s" % generation  : variance
    }
    if generation == 0:
        with open(log_file, "w") as file:
            yaml.dump(data, file, default_flow_style=False)
    else:
        with open(log_file, "a") as file:
            yaml.dump(data, file, default_flow_style=False)
