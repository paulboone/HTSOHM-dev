from time import sleep
import os

import yaml
from sqlalchemy import func
import sjs

sjs.load(os.path.join("settings","sjs.yaml"))
queue = sjs.get_job_queue()

from htsohm.utilities import write_config_file, read_config_file, evaluate_convergence
from htsohm.utilities import save_convergence
from htsohm.htsohm import seed_generation, queue_all_materials, queue_create_next_generation
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

def generation_write_complete(run_id, generation):
    materials_per_generation = read_config_file(run_id)['children-per-generation']
    materials_successfully_written =  session \
        .query(func.count(Material.id)) \
        .filter(
            run_id == run_id, Material.generation == generation,
            Material.write_check == 'done'
        ) \
        .all()[0][0]
    return materials_successfully_written == materials_per_generation

def manage_run(run_id, generation):
    config = read_config_file(run_id)
    if generation > config['max-number-of-generations']:
        print("Max generations exceeded; terminating run.")
        final_convergence = evaluate_convergence(run_id)
        save_convergence(run_id, generation, final_convergence)
        return -1 # -1 means we're done
    elif generation == 0:
        seed_generation(run_id, config['children-per-generation'],
            config['number-of-atom-types'])
        queue_all_materials(run_id, generation, queue)
        generation += 1
    elif generation >= 1:
        if not generation_write_complete(run_id, generation):
            convergence = evaluate_convergence(run_id)
            save_convergence(run_id, generation - 1, convergence)
            if convergence <= config['acceptance-value']:
                print('Desired convergence attained; terminating run.')
                return -1
            update_strength_array(run_id, generation)
            queue_create_next_generation(run_id, generation, queue)
        else:
            queue_all_materials(run_id, generation, queue)
            generation += 1
    return generation
