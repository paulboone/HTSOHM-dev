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
from htsohm.dummy_test import dummy_test
from htsohm.utilities import write_config_file, evaluate_convergence, save_convergence

def simulate_all_materials(run_id, generation):
    """simulate methane loading, helium void fraction, and surface area for seed population"""
    materials = session.query(Material).filter(Material.run_id == run_id, Material.generation == generation).all()
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
    materials = session.query(Material).filter(Material.run_id == run_id, Material.generation == generation).all()
    for material in materials:
        queue.enqueue(hpc_job_run_all_simulations, material.id, timeout=60*60)

def seed_generation(run_id, children_per_generation, number_of_atomtypes, queue=None):
    generation = 0
    write_seed_definition_files(run_id, children_per_generation, number_of_atomtypes)
    session.commit()
    if queue is not None:
        queue_all_materials(run_id, generation, queue)
    else:
        simulate_all_materials(run_id, generation)

def next_generation(run_id, children_per_generation, generation, queue=None):
    if generation == 1:
        create_strength_array(run_id)
    elif generation >= 2:
        recalculate_strength_array(run_id, generation)
    for i in range(children_per_generation):
        test_complete = False
        while not test_complete:
            parent_id = select_parent(run_id)
            test_complete = dummy_test(run_id, parent_id)
            session.commit()
        write_child_definition_files(run_id, generation, parent_id)
        session.commit()
    if queue is not None:
        queue_all_materials(run_id, generation, queue)
    else:
        simulate_all_materials(run_id, generation)

def htsohm(children_per_generation,    # number of materials per generation
           number_of_atomtypes,        # number of atom-types per material
           strength_0,                 # intial strength parameter
           number_of_bins,             # number of bins for analysis
           max_generations=20,         # maximum number of generations
           acceptance_value=-0.5):      # desired degree of `convergence`
    ############################################################################
    # write run-configuration file
    run_id = write_config_file(children_per_generation, number_of_atomtypes, strength_0,
        number_of_bins, max_generations)["run-id"]

    convergence = acceptance_value + 1          # initialize convergence with arbitrary value
    for generation in range(max_generations):
        while convergence >= acceptance_value:
            if generation == 0:                     # SEED GENERATION
                seed_generation(run_id, children_per_generation, number_of_atomtypes)
                convergence = evaluate_convergence(run_id)
                save_convergence(run_id, generation, convergence)
            elif generation >= 1:                   # FIRST GENERATION, AND ON...
                next_generation(run_id, children_per_generation, generation)
                convergence = evaluate_convergence(run_id)
                save_convergence(run_id, generation, convergence)
            print('convergence:\t%s' % convergence)
            generation += 1
