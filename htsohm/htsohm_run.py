
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


def matrix_of_empty_lists(shape):
    m = np.empty(shape, dtype=object)
    for i in np.ndindex(m.shape):
        m[i] = list([])
    return m

def load_restart_db(gen, num_bins, bin_ranges, config, session):
    mats = session.query(Material).filter(Material.generation < gen).all()
    ids = np.array([m.id for m in mats])
    props = np.array([property_tuple(m, config["bin_properties"]) for m in mats])

    bin_materials = matrix_of_empty_lists([num_bins]*len(config["bin_properties"]))

    bins = calc_bins(props, num_bins, bin_ranges)
    for i, bin_tuple in enumerate(bins):
        bin_materials[bin_tuple].append(i)

    return ids, props, bin_materials

def check_db_materials_for_restart(expected_num_materials, session, delete_excess=False):
    """Checks for if there are enough or extra materials in the database."""
    num_materials = session.query(Material).count()

    if num_materials < expected_num_materials:
        print("The database has fewer materials in it than the restart file indicated.")
        print("Is this the right database and restart file?")
        sys.exit(1)

    if num_materials > expected_num_materials:
        print("The database has an extra %d materials in it." % (num_materials - expected_num_materials))
        if (delete_excess):
            print("deleting from materials where id > %d" % expected_num_materials)
            db.delete_extra_materials(expected_num_materials)
        else:
            print("Is this the right database and restart file?")
            raise(Exception("db"))

    max_id = session.query(Material).order_by(Material.id.desc()).limit(1)[0].id
    if max_id != expected_num_materials:
        print("The max id of the materials table is not the same as the number of mateirals. Do "
            "you have missing rows?")
        sys.exit(1)

    return True

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

def property_tuple(material, properties):
    results = []
    for propname in properties:
        p = getattr(material, propname)
        if "log10" in properties[propname] and properties[propname]["log10"]:
            # minimum value of property is start of range
            min_p = 10.0 ** properties[propname]['range'][0]
            results += [math.log10(max(p, min_p))]
        else:
            results += [p]
    return tuple(results)

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
    return (material.id, property_tuple(material, config["bin_properties"]))


def parallel_simulate_generation(generator, num_processes, parent_ids, config, gen, children_per_generation):
    worker_metadata = (generator, config, gen)

    if parent_ids is None:
        parent_ids = [0] * (children_per_generation) # should only be needed for random!

    with Pool(processes=num_processes, initializer=init_worker, initargs=[worker_metadata]) as pool:
        results = pool.map(simulate_generation_worker, parent_ids)

    ids, props = zip(*results)
    return (np.array(ids), np.array(props))

def select_parents(children_per_generation, ids, props, bin_materials, config):
    if config['generator_type'] == 'random':
        return (None, [])
    elif config['selector_type'] == 'simplices-or-hull':
        return selector_tri.choose_parents(children_per_generation, ids, props, config['simplices_or_hull'])
    elif config['selector_type'] == 'density-bin':
        return selector_bin.choose_parents(children_per_generation, ids, props, bin_materials)
    elif config['selector_type'] == 'neighbor-bin':
        return selector_neighbor_bin.choose_parents(children_per_generation, ids, props, bin_materials)
    elif config['selector_type'] == 'best':
        return selector_best.choose_parents(children_per_generation, ids, props)
    elif config['selector_type'] == 'specific':
        return selector_specific.choose_parents(children_per_generation, ids, props, config['selector_specific_id'])


def htsohm_run(config_path, restart_generation=-1, override_db_errors=False, num_processes=1, max_generations=None, return_run_vars=False):
    """

    if return_run_vars, htsohm_run returns the main run vars, which are:
    - ids: the database ids. Since this is sqlite which is 1-indexed, this should always be one greater than the index.
    - props: a tuple of properties to bin per material (0 indexed)
    - bin_materials: an n-dimensional matrix (one dimension per bin dimension); typically accessed
        bin_materials[bin_tuple], where each entry is a list of material indices in the bin (0 indexed).
    """

    def _update_bin_materials(all_bins, start_index):
        nonlocal bin_materials
        for i, bin_tuple in enumerate(all_bins):
            bin_materials[bin_tuple].append(i + start_index)

    config = load_config_file(config_path)
    os.makedirs(config['output_dir'], exist_ok=True)
    print(config)

    children_per_generation = config['children_per_generation']

    bin_ranges = [config['bin_properties'][p]["range"] for p in config['bin_properties']]
    num_bins = config["num_bins"]
    if max_generations is None:
        max_generations = config['max_generations']

    engine, session = db.init_database(config["database_connection_string"], config['properties'],
                backup=(restart_generation > 1))

    print('{:%Y-%m-%d %H:%M:%S}'.format(datetime.now()))
    if restart_generation > 1:
        print("Restarting from database at generation: %s" % restart_generation)
        check_db_materials_for_restart((restart_generation - 1)*children_per_generation, session, delete_excess=override_db_errors)
        ids, props, bin_materials = load_restart_db(restart_generation, num_bins, bin_ranges, config, session)
        print("We loaded %d materials from the database" % len(props))
        start_gen = restart_generation
    else:
        if session.query(Material).count() > 0:
            print("ERROR: cannot have existing materials in the database for a new run")
            sys.exit(1)

        # generate initial generation of random materials
        print("applying random seed to initial points: %d" % config['initial_points_random_seed'])
        random.seed(config['initial_points_random_seed'])
        ids, props = parallel_simulate_generation(generator.random.new_material, num_processes, None,
                        config, gen=1, children_per_generation=config['children_per_generation'])
        random.seed() # flush the seed so that only the initial points are set, not generated points

        # setup initial bins
        bin_materials =  matrix_of_empty_lists([num_bins]*len(config["bin_properties"]))
        all_bins = calc_bins(props, num_bins, bin_ranges)
        _update_bin_materials(all_bins, 0)

        print([print(i, b, all_bins[i]) for i, b in enumerate(props)])
        start_gen = 2

    if config['generator_type'] == 'random':
        generator_method = generator.random.new_material
    elif config['generator_type'] == 'mutate':
        generator_method = generator.mutate.mutate_material

    # first generation is randomly determined above, so start_gen should always be >= 2.
    for gen in range(start_gen, max_generations + 1):
        # mutate materials and simulate properties
        parents_ids, parents_props = select_parents(children_per_generation, ids, props, bin_materials, config)
        new_ids, new_props = parallel_simulate_generation(generator_method, num_processes, parents_ids,
                                config, gen=gen, children_per_generation=config['children_per_generation'])

        # track bins
        all_bins = calc_bins(new_props, num_bins, bin_ranges)
        _update_bin_materials(all_bins, (gen - 1) * children_per_generation)

        # evaluate algorithm effectiveness
        explored_bins = [i for i,v in np.ndenumerate(bin_materials) if len(v) > 0]
        bin_fraction_explored = len(explored_bins) / num_bins ** 2
        print_block('GENERATION %s: %5.2f%%' % (gen, bin_fraction_explored * 100))

        ids = np.append(ids, new_ids, axis=0)
        props = np.append(props, new_props, axis=0)

    if return_run_vars:
        return ids, props, bin_materials
    return None

    # with open("pm.csv", 'w', newline='') as f:
    #     output_csv_from_db(session, output_file=f)
    #
    # with open("pm-binned.csv", 'w', newline='') as f:
    #     # column 12 is void_fraction_geo, 13 is methane loading
    #     csv_add_bin_column("pm.csv", [(12, *prop1range, num_bins), (13, *prop2range, num_bins)], output_file=f)
