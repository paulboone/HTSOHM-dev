# standard imports
import os

# related third party imports
import yaml

# local application/library specific imports

def load_config_file(file_name):
    """Reads input file.

    Input files must be in .yaml format, see htsohm.sample.yaml
    """

    with open(file_name) as file:
         config = yaml.load(file)
    return config

def save_convergence_file(run_id, generation, variance):
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
