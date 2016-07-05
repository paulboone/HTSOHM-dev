from time import sleep
import os

import yaml
from sqlalchemy import func
import sjs

sjs.load(os.path.join("settings","sjs.yaml"))
job_queue = sjs.get_job_queue()

from htsohm.utilities import write_config_file, read_config_file
from htsohm.htsohm import seed_generation, queue_all_materials, queue_create_next_gen
from htsohm.htsohm import update_strength_array
from htsohm.runDB_declarative import Material, session

def start_run(
        children_per_generation,    # number of materials per generation
        number_of_atomtypes,        # number of atom-types per material
        strength_0,                 # intial strength parameter
        number_of_bins,             # number of bins for analysis
        max_generations=100,        # maximum number of generations
        dummy_test_trials=3,        # number of re-simulations for dummy-test
        acceptance_value=-0.001):   # desired degree of `convergence`

    run_id = write_config_file(children_per_generation, number_of_atomtypes, strength_0,
        number_of_bins, max_generations, dummy_test_trials, acceptance_value)['run-id']

    return run_id

def manage_run(run_id, generation):
    config = read_config_file(run_id)
    
    # prepare, mutate, and queue up the next generation
    if generation > config['max-number-of-generations']:
        print("Max generations exceeded; terminating run.")
        return -1 # -1 means we're done
    elif generation == 0:
        # SEED GENERATION
        seed_generation(run_id,
                        config['children-per-generation'], config['number-of-atom-types'])
        queue_all_materials(run_id, generation, queue)
    elif generation >= 1:
        # FIRST GENERATION, AND ON...
        update_strength_array(run_id, generation)
        queue_create_next_gen(run_id, generation, queue)
        queue_all_materials(run_id, generation, queue)
    return generation + 1
