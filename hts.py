#!/usr/bin/env python3

import click
from sqlalchemy.sql import func

from htsohm.binning import select_parent
from htsohm.generate import write_seed_definition_files
from htsohm.mutate import write_child_definition_files, create_strength_array
from htsohm.db import Material, session
from htsohm.simulate import run_all_simulations
from htsohm.utilities import read_config_file, write_config_file

def materials_in_generation(run_id, generation):
    return session.query(Material).filter(
        Material.run_id == run_id,
        Material.generation == generation
    ).count()

def last_generation(run_id):
    return session.query(func.max(Material.generation)).filter(
        Material.run_id == run_id,
    )[0][0]

@click.group()
def hts():
    pass

@hts.command()
@click.argument('num_atomtypes', type=click.INT)
@click.argument('strength', type=click.FLOAT)
@click.argument('num_bins', type=click.INT)
@click.argument('children_in_generation', type=click.INT)
@click.argument('num_seeds', type=click.INT)
@click.argument('acceptance_value', type=click.FLOAT)
def start(num_atomtypes, strength, num_bins, children_in_generation, num_seeds, acceptance_value):
    config = write_config_file(num_atomtypes, strength, num_bins, children_in_generation, num_seeds, acceptance_value)
    run_id = config["run-id"]
    create_strength_array(run_id)
    print("Run created with id: %s" % run_id)

@hts.command()
@click.argument("run_id")
def launch_worker(run_id):
    config = read_config_file(run_id)

    gen = last_generation(run_id) or 0

    converged = False
    while not converged:
        size_of_generation = config['children-in-generation']
        if gen == 0 and config['num-seeds']:
            size_of_generation = config['num-seeds']

        while materials_in_generation(run_id, gen) < size_of_generation:
            if gen == 0:
                print("writing new seed...")
                material = write_seed_definition_files(run_id, config['number-of-atom-types'])
            else:
                print("creating / simulating new material")
                parent_id = select_parent(run_id, max_generation=(gen - 1),
                                                  generation_limit=config['children-in-generation'])
                material = write_child_definition_files(run_id, parent_id, gen)

            run_all_simulations(material)
            session.add(material)
            session.commit()

            material.generation_index = material.calculate_generation_index()
            if material.generation_index < config['children-in-generation']:
                session.add(material)
            else:
                # delete excess rows
                session.delete(material)
            session.commit()
        gen += 1

        # no convergance test at present!

if __name__ == '__main__':
    hts()
