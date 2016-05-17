import sys
import os
from datetime import datetime
import yaml

import numpy as np

from htsohm import generate as gen
from htsohm import simulate as sim
from htsohm import binning as bng
from htsohm import mutate as mut
from htsohm.runDB_declarative import RunData, session

def write_config_file(children_per_generation, number_of_atomtypes, strength_0,
    number_of_bins, max_generations):
    run_id = datetime.now().isoformat()
    wd = os.environ['HTSOHM_DIR']      # specify working directory
    config_file = os.path.join(wd, 'config', run_id + '.yaml')
    run_config = {
        "run-id" : run_id,
        "children-per-generation" : children_per_generation,
        "number-of-atom-types" : number_of_atomtypes,
        "initial-mutation-strength" : strength_0,
        "number-of-bins" : number_of_bins,
        "max-number-of-generations" : max_generations}
    with open(config_file, "w") as file:
        yaml.dump(run_config, file, default_flow_style=False)
    return run_config["run-id"]

def seed_generation(run_id, children_per_generation, number_of_atomtypes):
    generation = 0
    for material in range(children_per_generation):
        new_material = RunData(run_id, generation, 'none')
        session.add(new_material)
    session.commit()
    gen.write_seed_definition_files(run_id, children_per_generation, number_of_atomtypes)
    seed = session.query(RunData).filter(RunData.run_id == run_id, RunData.generation == 0).all()
    for material in seed:
        sim.run_all_simulations(material.id)
    session.commit()

def next_generation(run_id, children_per_generation, generation):
    s = session    # rename database objects to shorted queries
    db = RunData

    for material in range(children_per_generation):
        new_material = db(run_id, generation, 'none')
        s.add(new_material)
    s.commit()
    status = "Dummy test:   RUNNING"
    while status == "Dummy test:   RUNNING":
        # Select parents, add IDs to database...
        next_generation_list = bng.select_parents(run_id, children_per_generation, generation)
        s.commit()             # parent_ids staged by bng.select_parents()
        status = sim.dummy_test(run_id, next_generation_list, status, generation)
        s.commit()             # within loop so "fail" parents aren't re-selected
    if generation == 1:
        mut.create_strength_array(run_id)                  # Create strength-parameter array `run_id`.npy
    elif generation >= 2:
        mut.recalculate_strength_array(run_id, generation) # Recalculate strength-parameters, as needed
    mut.write_mutant_definition_files(run_id, generation)  # Create first generation of child-materials
    children = s.query(db).filter(db.run_id == run_id, db.generation == generation).all()
    for material in children:
        sim.run_all_simulations(material.id)
    s.commit()

def htsohm(children_per_generation,    # number of materials per generation
           number_of_atomtypes,        # number of atom-types per material
           strength_0,                 # intial strength parameter
           number_of_bins,             # number of bins for analysis
           max_generations=20):        # maximum number of generations

    run_id = write_config_file(children_per_generation, number_of_atomtypes, strength_0, 
        number_of_bins, max_generations)

    for generation in range(max_generations):
        if generation == 0:                     # SEED GENERATION
            seed_generation(run_id, children_per_generation, number_of_atomtypes)
        elif generation >= 1:                   # FIRST GENERATION, AND ON...
            next_generation(run_id, children_per_generation, generation)
