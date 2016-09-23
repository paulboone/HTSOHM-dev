import os

from htsohm.files import load_config_file

config = {}

def _init(run_id):
    wd = os.environ['HTSOHM_DIR']
    config_file = os.path.join(wd, 'config', run_id + '.yaml')
    config.update(load_config_file(config_file))
    return config
