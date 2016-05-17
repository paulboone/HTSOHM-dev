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
    """Writes run-specific parameters to /config/<run_id>.yaml.

    This function writes run-specific parameters to a configuration file that is loaded/written
    at different stages throughout the overall HTSOHM routine.
    """
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

    ############################################################################
    # initialize materials in database
    for material in range(children_per_generation):
        new_material = RunData(run_id, generation, 'none')
        session.add(new_material)
    session.commit()

    ############################################################################
    # write material-definition files for seed population:
    #     ___________________________________________________________________________________
    #       $MAT_DIR/<material-name>.cif                          |  structural information
    #     ________________________________________________________|__________________________
    #       $FF_DIR/<material-name>/force_field.def               |  overwrite LJ-parameters
    #       $FF_DIR/<material-name>/force_field_mixing_rules.def  |  LJ-parameters
    #       $FF_DIR/<material-name>/pseudo_atoms.def              |  atomic charges, masses
    #     ________________________________________________________|__________________________
    gen.write_seed_definition_files(run_id, children_per_generation, number_of_atomtypes)

    ############################################################################
    # simulate methane loading, helium void fraction, and surface area for seed population
    seed = session.query(RunData).filter(RunData.run_id == run_id, RunData.generation == 0).all()
    for material in seed:
        sim.run_all_simulations(material.id)
    session.commit()

def next_generation(run_id, children_per_generation, generation):
    s = session    # rename database objects to shorted queries
    db = RunData

    ############################################################################
    # initialize materials in database
    for material in range(children_per_generation):
        new_material = db(run_id, generation, 'none')
        s.add(new_material)
    s.commit()

    ############################################################################
    # select potential parent-materials and run them in dummy-test
    status = "Dummy test:   RUNNING"
    while status == "Dummy test:   RUNNING":
        next_generation_list = bng.select_parents(run_id, children_per_generation, generation)
        s.commit()             # parent_ids added to database
        status = sim.dummy_test(run_id, next_generation_list, status, generation)
        s.commit()             # dummy_test_results added to database

    ############################################################################
    # selected parent-materials are mutated to create next generation
    if generation == 1:
        mut.create_strength_array(run_id)                  # create strength-parameter array `run_id`.npy
    elif generation >= 2:
        mut.recalculate_strength_array(run_id, generation) # recalculate strength-parameters, as needed
    mut.write_mutant_definition_files(run_id, generation)  # create child-materials

    ############################################################################
    # simulate methane loading, helium void fraction, and surface area for seed population
    children = s.query(db).filter(db.run_id == run_id, db.generation == generation).all()
    for material in children:
        sim.run_all_simulations(material.id)
    s.commit()

def htsohm(children_per_generation,    # number of materials per generation
           number_of_atomtypes,        # number of atom-types per material
           strength_0,                 # intial strength parameter
           number_of_bins,             # number of bins for analysis
           max_generations=20):        # maximum number of generations

    ############################################################################
    # write run-configuration file
    run_id = write_config_file(children_per_generation, number_of_atomtypes, strength_0,
        number_of_bins, max_generations)

    for generation in range(max_generations):
        if generation == 0:                     # SEED GENERATION
            seed_generation(run_id, children_per_generation, number_of_atomtypes)
        elif generation >= 1:                   # FIRST GENERATION, AND ON...
            next_generation(run_id, children_per_generation, generation)
