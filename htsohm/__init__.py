import os

import htsohm
from htsohm.files import load_config_file

config = {}

def _init(run_id):
    htsohm_dir = os.path.dirname(os.path.dirname(htsohm.__file__))
    config_file = os.path.join(htsohm_dir, 'config', run_id + '.yaml')
    config.update(load_config_file(config_file))
    return config
