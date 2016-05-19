import sys
import os

import numpy as np

from htsohm.utilities import write_config_file
from htsohm import generate as gen
from htsohm import simulate as sim
from htsohm import binning as bng
from htsohm import mutate as mut
from htsohm.runDB_declarative import Material, session

def init_materials_in_database(run_id, children_per_generation, generation):
    """initialize materials in database with run_id and generation"""

    for material in range(children_per_generation):
        new_material = Material(run_id, generation, 'none')
        session.add(new_material)
    session.commit()

def simulate_all_materials(run_id, generation):
    """simulate methane loading, helium void fraction, and surface area for seed population"""
    materials = session.query(Material).filter(Material.run_id == run_id, Material.generation == generation).all()
    for material in materials:
        sim.run_all_simulations(material.id)
    session.commit()

def hpc_job_run_all_simulations(material_id):
    sim.run_all_simulations(material_id)
    session.commit()

def queue_all_materials(run_id, generation, queue):
    """same as simulate_all_materials, except queues the jobs in the job server"""
    materials = session.query(Material).filter(Material.run_id == run_id, Material.generation == generation).all()
    for material in materials:
        queue.enqueue(hpc_job_run_all_simulations, material.id, timeout=60*60)

def screen_parents(run_id, children_per_generation, generation):
    """select potential parent-materials and run them in dummy-test"""
    test_complete = False
    while not test_complete:
        next_generation_list = bng.select_parents(run_id, children_per_generation, generation)
        session.commit()             # parent_ids added to database
        test_complete = sim.dummy_test(run_id, next_generation_list, generation)
        session.commit()             # dummy_test_results added to database

def create_next_generation(run_id, generation):
    """once screened, parent-materials are mutated to create next generation"""
    if generation == 1:
        mut.create_strength_array(run_id)                  # create strength-parameter array `run_id`.npy
    elif generation >= 2:
        mut.recalculate_strength_array(run_id, generation) # recalculate strength-parameters, as needed
    mut.write_children_definition_files(run_id, generation)      # create child-materials

def seed_generation(run_id, children_per_generation, number_of_atomtypes, queue=None):
    generation = 0
    init_materials_in_database(run_id, children_per_generation, generation)
    gen.write_seed_definition_files(run_id, children_per_generation, number_of_atomtypes)
    if queue is not None:
        queue_all_materials(run_id, generation, queue)
    else:
        simulate_all_materials(run_id, generation)

def next_generation(run_id, children_per_generation, generation, queue=None):
    init_materials_in_database(run_id, children_per_generation, generation)
    screen_parents(run_id, children_per_generation, generation)
    create_next_generation(run_id, generation)
    if queue is not None:
        queue_all_materials(run_id, generation, queue)
    else:
        simulate_all_materials(run_id, generation)

def htsohm(children_per_generation,    # number of materials per generation
           number_of_atomtypes,        # number of atom-types per material
           strength_0,                 # intial strength parameter
           number_of_bins,             # number of bins for analysis
           max_generations=20):        # maximum number of generations

    ############################################################################
    # write run-configuration file
    run_id = write_config_file(children_per_generation, number_of_atomtypes, strength_0,
        number_of_bins, max_generations)["run-id"]

    seed_generation(run_id, children_per_generation, number_of_atomtypes)
    for generation in range(1,max_generations):
        next_generation(run_id, children_per_generation, generation)
