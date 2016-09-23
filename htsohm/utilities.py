# stanard imports
import os
from math import sqrt

# related third party imports
from datetime import datetime
import yaml
from sqlalchemy import func

# local application/library specific imports
from htsohm.db import Base, Material, session

def load_config_file(file_name):
    """Reads input file.

    Input files must be in .yaml format, see input_file.sample.yaml
    """

    with open(file_name) as file:
         config = yaml.load(file)
    return config


def evaluate_convergence(run_id):
    '''Counts number of materials in each bin and returns variance of these counts.'''
    bin_counts = session \
        .query(func.count(Material.id)) \
        .filter(Material.run_id == run_id) \
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
