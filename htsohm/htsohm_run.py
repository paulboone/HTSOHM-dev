
from datetime import datetime
from glob import glob
import math
from multiprocessing import Pool
import os
import random
import shutil
import sys

import numpy as np
from sqlalchemy.orm import joinedload

from htsohm import generator, load_config_file, db
from htsohm.bins import calc_bins
from htsohm.bin.output_csv import output_csv_from_db, csv_add_bin_column
from htsohm.db import Material
from htsohm.simulation.run_all import run_all_simulations
import htsohm.select.triangulation as selector_tri
import htsohm.select.density_bin as selector_bin
import htsohm.select.best as selector_best
import htsohm.select.specific as selector_specific
import htsohm.select.neighbor_bin as selector_neighbor_bin
from htsohm.slog import init_slog, get_slog, slog

def print_block(string):
    print('{0}\n{1}\n{0}'.format('=' * 80, string))


def empty_lists_2d(x,y):
    return [[[] for j in range(x)] for i in range(y)]

def load_restart_db(gen, num_bins, prop1range, prop2range, session):
    mats = session.query(Material).options(joinedload("void_fraction"), joinedload("gas_loading")) \
                    .filter(Material.generation <= gen).all()
    box_d = np.array([m.id for m in mats])
    box_r = np.array([(m.void_fraction[0].get_void_fraction(),
                        math.log10(max(m.henrys_coefficient[0].co2_henrysv, 10.0 ** prop2range[0]))) for m in mats])

    bin_counts = np.zeros((num_bins, num_bins))
    bin_materials = empty_lists_2d(num_bins, num_bins)

    bins = calc_bins(box_r, num_bins, prop1range=prop1range, prop2range=prop2range)
    for i, (bx, by) in enumerate(bins):
        bin_counts[bx,by] += 1
        bin_materials[bx][by].append(i)

    start_gen = gen + 1
    return box_d, box_r, bin_counts, bin_materials, set(bins), start_gen

def check_db_materials_for_restart(expected_num_materials, session, delete_excess=False):
    """Checks for if there are enough or extra materials in the database."""
    extra_materials = session.query(Material).filter(Material.id > expected_num_materials).all()
    if len(extra_materials) > 0:
        print("The database has an extra %d materials in it." % len(extra_materials))
        if (delete_excess):
            print("deleting from materials where id > %d" % expected_num_materials)
            db.delete_extra_materials(expected_num_materials)
        else:
            print("Is this the right database and restart file?")
            sys.exit(1)

    num_materials = session.query(Material).count()
    if num_materials < expected_num_materials:
        print("The database has fewer materials in it than the restart file indicated.")
        print("Is this the right database and restart file?")
        sys.exit(1)

def init_worker(worker_metadata):
    """initialization function for worker that inits the database and gets a worker-specific
    session."""
    global config
    global generator
    global gen
    global worker_session
    generator, config, gen = worker_metadata
    _, worker_session = db.init_database(config["database_connection_string"], config["properties"])
    return

def simulate_generation_worker(parent_id):
    """gets most of its parameters from the global worker_metadata set in the
    parallel_simulate_generation method."""
    init_slog()
    if parent_id > 0:
        parent = worker_session.query(Material).get(int(parent_id))
        material = generator(parent, config["structure_parameters"])
    else:
        material = generator(config["structure_parameters"])

    run_all_simulations(material, config)
    material.generation = gen
    worker_session.add(material)
    worker_session.commit()

    print(get_slog())
    return (material.id, (material.void_fraction[0].get_void_fraction(),
                          math.log10(max(material.henrys_coefficient[0].co2_henrysv, 10.0 ** config['prop2range'][0]))))

def parallel_simulate_generation(generator, num_processes, parent_ids, config, gen, children_per_generation):
    worker_metadata = (generator, config, gen)

    if parent_ids is None:
        parent_ids = [0] * (children_per_generation) # should only be needed for random!

    with Pool(processes=num_processes, initializer=init_worker, initargs=[worker_metadata]) as pool:
        results = pool.map(simulate_generation_worker, parent_ids)

    box_d, box_r = zip(*results)
    return (np.array(box_d), np.array(box_r))

def select_parents(children_per_generation, box_d, box_r, bin_materials, config):
    if config['generator_type'] == 'random':
        return (None, [])
    elif config['selector_type'] == 'simplices-or-hull':
        return selector_tri.choose_parents(children_per_generation, box_d, box_r, config['simplices_or_hull'])
    elif config['selector_type'] == 'density-bin':
        return selector_bin.choose_parents(children_per_generation, box_d, box_r, bin_materials)
    elif config['selector_type'] == 'neighbor-bin':
        return selector_neighbor_bin.choose_parents(children_per_generation, box_d, box_r, bin_materials)
    elif config['selector_type'] == 'best':
        return selector_best.choose_parents(children_per_generation, box_d, box_r)
    elif config['selector_type'] == 'specific':
        return selector_specific.choose_parents(children_per_generation, box_d, box_r, config['selector_specific_id'])


def htsohm_run(config_path, restart_generation=-1, override_db_errors=False, num_processes=1, max_generations=None):

    def _update_bins_counts_materials(all_bins, bins, start_index):
        nonlocal bin_counts, bin_materials
        for i, (bx, by) in enumerate(all_bins):
            bin_counts[bx,by] += 1
            bin_materials[bx][by].append(i + start_index)
        new_bins = set(all_bins) - bins
        return new_bins, bins.union(new_bins)


    config = load_config_file(config_path)
    os.makedirs(config['output_dir'], exist_ok=True)
    print(config)

    children_per_generation = config['children_per_generation']
    prop1range = config['prop1range']
    prop2range = config['prop2range']
    num_bins = config['num_bins']

    if max_generations is None:
        max_generations = config['max_generations']

    engine, session = db.init_database(config["database_connection_string"], config['properties'],
                backup=(restart_generation > 0))

    print('{:%Y-%m-%d %H:%M:%S}'.format(datetime.now()))
    if restart_generation >= 0:
        print("Restarting from database using generation: %s" % restart_generation)
        box_d, box_r, bin_counts, bin_materials, bins, start_gen = load_restart_db(
            restart_generation, num_bins, prop1range, prop2range, session)

        print("Restarting at generation %d\nThere are currently %d materials" % (start_gen, len(box_r)))
        check_db_materials_for_restart(len(box_r), session, delete_excess=override_db_errors)
    else:
        if session.query(Material).count() > 0:
            print("ERROR: cannot have existing materials in the database for a new run")
            sys.exit(1)

        # generate initial generation of random materials
        print("applying random seed to initial points: %d" % config['initial_points_random_seed'])
        random.seed(config['initial_points_random_seed'])
        box_d, box_r = parallel_simulate_generation(generator.random.new_material, num_processes, None,
                        config, gen=0, children_per_generation=config['children_per_generation'])
        random.seed() # flush the seed so that only the initial points are set, not generated points

        # setup initial bins
        bin_counts = np.zeros((num_bins, num_bins))
        bin_materials = empty_lists_2d(num_bins, num_bins)
        all_bins = calc_bins(box_r, num_bins, prop1range=prop1range, prop2range=prop2range)
        new_bins, bins = _update_bins_counts_materials(all_bins, set(), 0)

        print([print(i, b, all_bins[i]) for i, b in enumerate(box_r)])
        start_gen = 1

    if config['generator_type'] == 'random':
        generator_method = generator.random.new_material
    elif config['generator_type'] == 'mutate':
        generator_method = generator.mutate.mutate_material

    for gen in range(start_gen, max_generations + 1):
        # mutate materials and simulate properties
        parents_d, parents_r = select_parents(children_per_generation, box_d, box_r, bin_materials, config)
        new_box_d, new_box_r = parallel_simulate_generation(generator_method, num_processes, parents_d,
                                config, gen=gen, children_per_generation=config['children_per_generation'])

        # track bins
        all_bins = calc_bins(new_box_r, num_bins, prop1range=prop1range, prop2range=prop2range)
        [print(i, b, all_bins[i]) for i, b in enumerate(new_box_r)]
        new_bins, bins = _update_bins_counts_materials(all_bins, bins, gen * children_per_generation)

        # evaluate algorithm effectiveness
        bin_fraction_explored = len(bins) / num_bins ** 2
        print_block('GENERATION %s: %5.2f%%' % (gen, bin_fraction_explored * 100))

        box_d = np.append(box_d, new_box_d, axis=0)
        box_r = np.append(box_r, new_box_r, axis=0)

    # with open("pm.csv", 'w', newline='') as f:
    #     output_csv_from_db(session, output_file=f)
    #
    # with open("pm-binned.csv", 'w', newline='') as f:
    #     # column 12 is void_fraction_geo, 13 is methane loading
    #     csv_add_bin_column("pm.csv", [(12, *prop1range, num_bins), (13, *prop2range, num_bins)], output_file=f)
