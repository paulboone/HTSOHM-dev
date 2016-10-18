#!/usr/bin/env python3
from datetime import datetime
import os

import click
import yaml
import RASPA2

import htsohm
from htsohm.files import load_config_file
from htsohm.htsohm import worker_run_loop

@click.group()
def hts():
    pass

@hts.command()
@click.argument('config_path',type=click.Path())
def start(config_path):
    config = load_config_file(config_path)
    htsohm_dir = os.path.dirname(os.path.dirname(htsohm.__file__))
    run_id = datetime.now().isoformat()
    config['run_id'] = run_id
    config['raspa2_dir'] = os.path.dirname(RASPA2.__file__)
    config['htsohm_dir'] = htsohm_dir
    
    run_dir = os.path.join(htsohm_dir, run_id)
    os.makedirs(run_dir, exist_ok=True)
    config_file = os.path.join(run_dir, 'run_parameters.yaml')
    with open(config_file, 'w') as file:
        yaml.dump(config, file, default_flow_style=False)
    print('Run created with id: %s' % run_id)

@hts.command()
@click.argument('run_id')
def launch_worker(run_id):
    htsohm._init(run_id)
    worker_run_loop(run_id)

if __name__ == '__main__':
    hts()
