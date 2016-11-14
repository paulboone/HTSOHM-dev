import os

import htsohm
from htsohm.files import load_config_file

config = {}

def _init(run_id):
    """Load config for run.

    Args:
        run_id (str): identification string for run.

    Returns:
        config (dict): parameters specified in config.

    """
    htsohm_dir = os.path.dirname(os.path.dirname(htsohm.__file__))
    run_dir = os.path.join(htsohm_dir, run_id)
    config_file = os.path.join(run_dir, 'config.yaml')
    config.update(load_config_file(config_file))
    return config
