#!/usr/bin/env python3
from datetime import datetime
import os

import click
import yaml

import htsohm
from htsohm.utilities import load_config_file
from htsohm.htsohm import worker_run_loop

@click.group()
def hts():
    pass

@hts.command()
@click.argument("config_path",type=click.Path())
def start(config_path):
    config = load_config_file(config_path)

    run_id = datetime.now().isoformat()
    config['run-id'] = run_id

    wd = os.environ['HTSOHM_DIR']
    config_file = os.path.join(wd, 'config', run_id + '.yaml')
    with open(config_file, "w") as file:
        yaml.dump(config, file, default_flow_style=False)

    print("Run created with id: %s" % run_id)

@hts.command()
@click.argument("run_id")
def launch_worker(run_id):
    htsohm._init(run_id)
    worker_run_loop(run_id)

if __name__ == '__main__':
    hts()
