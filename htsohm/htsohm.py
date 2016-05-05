# /usr/bin/env python

import sys
import os
import datetime

import numpy as np

from htsohm import generate as gen
from htsohm import simulate as sim
from htsohm import binning as bng
from htsohm import mutate as mut


def write_run_file(children_per_generation,
                   number_of_atomtypes,
                   strength_0,
                   number_of_bins,
                   max_generations):

    sys.path.insert(0, os.environ['HTSOHM_DIR'] + '/htsohm')

    start = datetime.datetime.now()
    run_id = ( "%s.%s.%s_%s.%s.%s" %
               (start.day, start.month, start.year,
                start.hour, start.minute, start.second) )

    wd = os.environ['HTSOHM_DIR']      # specify working directory

    run_file = open( wd + '/' + run_id + '.txt', "w")
    run_file.write( "Date:\t\t\t\t%s:%s:%s\n" % (start.day, start.month,
                                                 start.year) +
                    "Time:\t\t\t\t%s:%s:%s\n" % (start.hour, start.minute,
                                                 start.second) +
                    "Children per generation:\t%s\n" % (
                                                 children_per_generation) +
                    "Number of atom-types:\t\t%s\n" % (number_of_atomtypes) +
                    "Initial mutation strength:\t%s\n" % (strength_0) +
                    "Number of bins:\t\t\t%s\n" % (number_of_bins))
    run_file.close()

    return run_id


def grep_run_file(run_id):

    wd = os.environ['HTSOHM_DIR']
    with open( wd + '/' + run_id + '.txt' ) as origin:
        for line in origin:
            if "Children per generation:" in line:
                children_per_generation = int( line.split()[3] )
            elif "Number of atom-types:" in line:
                number_of_atomtypes = int( line.split()[3] )

    return children_per_generation, number_of_atomtypes


def seed_generation(run_id):

    children_per_generation, number_of_atomtypes = grep_run_file(run_id)
    generation = 0
    ids = gen_ids(generation, children_per_generation)

    gen.generate(children_per_generation, number_of_atomtypes, run_id)
    sim.add_rows(run_id, generation, ids)
    sim.simulate(run_id, ids)


def first_generation(run_id, strength_0):

    children_per_generation, number_of_atomtypes = grep_run_file(run_id)
    generation = 1
    ids = gen_ids(generation, children_per_generation)

    status = "Dummy test:   RUNNING"
    while status == "Dummy test:   RUNNING":
        # Select parents, add IDs to database...
        next_materials_list = bng.select_parents(run_id,
                                                 children_per_generation,
                                                 generation)
        sim.add_rows(run_id, generation, ids)
        bng.add_parent_ids(run_id, next_materials_list)
        status = sim.dummy_test(run_id,
                                next_materials_list,
                                status,
                                generation)
    mut.first_s(run_id, strength_0)    # Create strength-parameter array `run_id`.npy
    mut.mutate(run_id, generation)     # Create first generation of child-materials
    sim.simulate( run_id, gen_ids(generation, children_per_generation) )


def next_generation(run_id, generation):

    children_per_generation, number_of_atomtypes = grep_run_file(run_id)
    ids = gen_ids(generation, children_per_generation)

    status = "Dummy test:   RUNNING"
    while status == "Dummy test:   RUNNING":
        # Select parents, add IDs to database...
        next_materials_list = bng.select_parents(run_id,
                                                 children_per_generation,
                                                 generation)
        sim.add_rows(run_id, generation, ids)
        bng.add_parent_ids(run_id, next_materials_list)
        status = sim.dummy_test(run_id,
                                next_materials_list,
                                status,
                                generation)

    mut.calculate_s(run_id, generation)
    mut.mutate(run_id, generation)
    sim.simulate( run_id, ids )


def htsohm(children_per_generation,    # number of materials per generation
           number_of_atomtypes,        # number of atom-types per material
           strength_0,                 # intial strength parameter
           number_of_bins,             # number of bins for analysis
           max_generations=20):        # maximum number of generations

    run_id = write_run_file(children_per_generation,
                            number_of_atomtypes,
                            strength_0,
                            number_of_bins,
                            max_generations)

    for i in range(max_generations):
        if i == 0:                     # SEED GENERATION
            seed_generation(run_id)
        elif i == 1:                   # FIRST GENERATION
            first_generation(run_id, strength_0)
        elif i >= 2:                   # SECOND GENERATION(S), and on...
            next_generation(run_id, i)


def gen_ids(generation, children_per_generation):

    first = generation * children_per_generation
    last = (generation + 1) * children_per_generation
    gen_ids = np.arange(first, last)

    return gen_ids
