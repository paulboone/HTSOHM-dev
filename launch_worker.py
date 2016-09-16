#!/usr/bin/env python3

import click

from htsohm.binning import select_parent
from htsohm.generate import write_seed_definition_files
from htsohm.mutate import write_child_definition_files, create_strength_array
from htsohm.runDB_declarative import Material, session
from htsohm.simulate import run_all_simulations
from htsohm.utilities import read_config_file, write_config_file

def seeds_in_run(run_id):
    return session.query(Material).filter(
        Material.run_id == run_id,
        Material.seed == True
    ).count()

@click.group()
def hts():
    pass

@hts.command()
@click.argument('num_atomtypes', type=click.INT)
@click.argument('strength', type=click.FLOAT)
@click.argument('num_bins', type=click.INT)
@click.argument('num_seeds', type=click.INT)
@click.argument('acceptance_value', type=click.FLOAT)
def start(num_atomtypes, strength, num_bins, num_seeds, acceptance_value):
    config = write_config_file(num_atomtypes, strength, num_bins, num_seeds, acceptance_value)
    run_id = config["run-id"]
    create_strength_array(run_id)
    print("Run created with id: %s" % run_id)

@hts.command()
@click.argument("run_id")
def launch_worker(run_id):
    config = read_config_file(run_id)

    while seeds_in_run(run_id) < config['num-seeds']:
        print("writing new seed...")
        material = write_seed_definition_files(run_id, config['number-of-atom-types'])
        run_all_simulations(material)
        session.add(material)
        session.commit()

    converged = False
    while not converged:
        print("creating / simulating new material")
        parent_id = select_parent(run_id)
        material = write_child_definition_files(run_id, parent_id)
        run_all_simulations(material)
        session.add(material)
        session.commit()

        # no convergance test at present!

if __name__ == '__main__':
    hts()
