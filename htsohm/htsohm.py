# standard library imports
import os
import sys

# related third party imports
import numpy as np

# local application/library specific imports
from htsohm.binning import select_parent
from htsohm.generate import write_seed_definition_files
from htsohm.mutate import create_strength_array, recalculate_strength_array
from htsohm.mutate import write_child_definition_files
from htsohm.runDB_declarative import Material, session
from htsohm.simulate import run_all_simulations
from htsohm.dummy_test import screen_parent
from htsohm.utilities import write_config_file, evaluate_convergence, save_convergence
from htsohm.utilities import read_config_file

def simulate_all_materials(run_id, generation):
    """simulate methane loading, helium void fraction, and surface area for seed population"""
    materials = session \
        .query(Material) \
        .filter(
            Material.run_id == run_id, Material.generation == generation,
            Material.write_check == 'done'
        ).all()
    for material in materials:
        run_all_simulations(material.id)
    session.commit()

def hpc_job_run_all_simulations(material_id):
    print("======================================================================================")
    print("== manhpc_job_run_all_simulations %s" % material_id)

    run_all_simulations(material_id)
    session.commit()
    print("======================================================================================")

def queue_all_materials(run_id, generation, queue):
    """same as simulate_all_materials, except queues the jobs in the job server"""
    materials = session \
        .query(Material) \
        .filter(
            Material.run_id == run_id, Material.generation == generation,
            Material.write_check == 'done'
        ).all()
    for material in materials:
        queue.enqueue(hpc_job_run_all_simulations, material.id, timeout=60*60)

def create_next_generation(run_id, generation):
    children_per_generation = read_config_file(run_id)['children-per-generation']
    for i in range(children_per_generation):
        parent_id = screen_parent(run_id)
        write_child_definition_files(run_id, generation, parent_id)
    session.commit()

def hpc_job_create_material(run_id, generation):
    print("======================================================================================")
    print("== manhpc_job_dummy_test %s" % material_id)

    parent_id = screen_parent(run_id)
    write_child_definition_files(run_id, generation, parent_id)
    session.commit()
    print("======================================================================================")

def queue_create_next_generation(run_id, generation, queue):
    """same as create_next_generation, except queues the jobs in the job server"""
    children_per_generation = read_config_file(run_id)['children-per-generation']
    for i in range(children_per_generation):
        trials = read_config_file(run_id)['dummy-test-trials']
        queue.enqueue(hpc_job_create_material, run_id, generation, timeout=trials*60*60)

def seed_generation(run_id, children_per_generation, number_of_atomtypes):
    generation = 0
    write_seed_definition_files(run_id, children_per_generation, number_of_atomtypes)
    session.commit()

def update_strength_array(run_id, generation):
    if generation == 1:
        create_strength_array(run_id)
    elif generation >= 2:
        recalculate_strength_array(run_id, generation)

def htsohm(children_per_generation,    # number of materials per generation
           number_of_atomtypes,        # number of atom-types per material
           strength_0,                 # intial strength parameter
           number_of_bins,             # number of bins for analysis
           max_generations=20,         # maximum number of generations
           dummy_test_trials=1,        # number of re-simulations for dummy-test
           acceptance_value=-0.005):      # desired degree of `convergence`
    run_id = write_config_file(children_per_generation, number_of_atomtypes, strength_0,
        number_of_bins, max_generations, dummy_test_trials, acceptance_value)["run-id"]
    for generation in range(max_generations):
            if generation == 0:                     # SEED GENERATION
                seed_generation(run_id, children_per_generation, number_of_atomtypes)
                simulate_all_materials(run_id, generation)
            elif generation >= 1:                   # FIRST GENERATION, AND ON...
                convergence = evaluate_convergence(run_id)
                save_convergence(run_id, generation - 1, convergence)
                print('convergence:\t%s' % convergence)
                if convergence <= acceptance_value:
                    print('Desired convergence attained; terminating run.')
                    break
                update_strength_array(run_id, generation)
                create_next_generation(run_id, generation)
                simulate_all_materials(run_id, generation)
            generation += 1
