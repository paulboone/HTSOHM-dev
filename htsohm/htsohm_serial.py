
import os
from datetime import datetime
import math
import random

import numpy as np

from htsohm import pseudomaterial_generator, load_config_file, db
from htsohm.db import Material
from htsohm.simulation.run_all import run_all_simulations

from htsohm.figures import delaunay_figure
import htsohm.select.triangulation as selector_tri
import htsohm.select.density_bin as selector_bin


def print_block(string):
    print('{0}\n{1}\n{0}'.format('=' * 80, string))


def calc_bin(value, bound_min, bound_max, bins):
    """Find bin in parameter range.
    Args:
        value (float): some value, the result of a simulation.
        bound_min (float): lower limit, defining the parameter-space.
        bound_max (float): upper limit, defining the parameter-space.
        bins (int): number of bins used to subdivide parameter-space.
    Returns:
        Bin(int) corresponding to the input-value.
    """
    step = (bound_max - bound_min) / bins
    assigned_bin = (value - bound_min) // step
    assigned_bin = min(assigned_bin, bins-1)
    assigned_bin = max(assigned_bin, 0)
    return int(assigned_bin)

def calc_bins(box_r, num_bins, prop1range=(0.0, 1.0), prop2range=(0.0, 1.0)):
    return [(calc_bin(b[0], *prop1range, num_bins), calc_bin(b[1], *prop2range, num_bins)) for b in box_r]

def empty_lists_2d(x,y):
    return [[[] for j in range(x)] for i in range(y)]

def serial_runloop(config_path):
    """
    Args:
        run_id (str): identification string for run.

    """

    run_id = datetime.now().isoformat()
    config = load_config_file(config_path)
    os.makedirs(config['output_dir'], exist_ok=True)
    db.init_database(config["database_connection_string"])
    session = db.get_session()

    num_bins = config['number_of_convergence_bins']
    bin_counts = np.zeros((num_bins, num_bins))
    bin_materials = empty_lists_2d(num_bins, num_bins)

    children_per_generation = config['children_per_generation']
    # prop1 = config['prop1']
    # prop2 = config['prop2']
    prop1range = config['prop1range']
    prop2range = config['prop2range']

    print('{:%Y-%m-%d %H:%M:%S}'.format(datetime.now()))

    verbose = config['verbose'] if 'verbose' in config else False
    benchmarks = config['benchmarks']
    next_benchmark = benchmarks.pop(0)

    print(config)

    # generate initial generation of random materials
    box_d = np.zeros(children_per_generation, dtype=int)
    box_r = -1 * np.ones((children_per_generation, 2))

    if config['initial_points_random_seed']:
        print("applying random seed to initial points: %d" % config['initial_points_random_seed'])
        random.seed(config['initial_points_random_seed'])

    for i in range(children_per_generation):
        print("Material Index: ", i)
        material = pseudomaterial_generator.random.new_material(run_id, config["structure_parameters"])
        run_all_simulations(material, config)
        session.add(material)
        session.commit()

        box_d[i] = material.id
        box_r[i,:] = (material.void_fraction[0].void_fraction, material.gas_loading[0].absolute_volumetric_loading)
        # box_r[i,:] = (material[prop1], material[prop2])

    random.seed() # flush the seed so that only the initial points are set, not generated points

    all_bins = calc_bins(box_r, num_bins, prop1range=prop1range, prop2range=prop2range)
    for i, (bx, by) in enumerate(all_bins):
        bin_counts[bx,by] += 1
        bin_materials[bx][by].append(i)
    bins = set(all_bins)
    print("bins", bins)

    output_path = os.path.join(config['output_dir'], "triplot_0.png")
    delaunay_figure(box_r, num_bins, output_path, bins=bin_counts, \
                    title="Starting random materials", \
                    prop1range=prop1range, prop2range=prop2range)


    last_benchmark_reached = False

    for gen in range(1, config['max_generations'] + 1):

        parents_r = parents_d = []
        perturbation_methods = None
        if config['selector_type'] == 'simplices-or-hull':
            parents_d, parents_r = selector_tri.choose_parents(children_per_generation, box_d, box_r, config['simplices_or_hull'])
            perturbation_methods = [""] * children_per_generation
        elif config['selector_type'] == 'density-bin':
            parents_d, parents_r = selector_bin.choose_parents(children_per_generation, box_d, box_r, bin_materials)
            print("parents_d: ", parents_d)
            perturbation_methods = [""] * children_per_generation

        # mutate materials and simulate properties
        new_box_d = np.zeros(children_per_generation)
        new_box_r = -1 * np.ones((children_per_generation, 2))
        for i in range(children_per_generation):
            print("Material Index: ", i + gen * children_per_generation)
            if config['generator_type'] == 'random':
                material = pseudomaterial_generator.random.new_material(run_id, config["structure_parameters"])
            elif config['generator_type'] == 'mutate':
                print("parents_d[i]: ", parents_d[i])
                material = pseudomaterial_generator.mutate.mutate_material(run_id, parents_d[i], config["structure_parameters"])
                perturbation_methods[i] = material.perturbation

            run_all_simulations(material, config)
            session.add(material)
            session.commit()

            new_box_d[i] = material.id
            new_box_r[i,:] = (material.void_fraction[0].void_fraction, material.gas_loading[0].absolute_volumetric_loading)
            # new_box_r[i,:] = (material[prop1], material[prop2])

        # TODO: bins for methane loading?
        all_bins = calc_bins(new_box_r, num_bins, prop1range=prop1range, prop2range=prop2range)
        for i, (bx, by) in enumerate(all_bins):
            bin_counts[bx,by] += 1
            material_index = i + gen * children_per_generation
            bin_materials[bx][by].append(material_index)
        new_bins = set(all_bins) - bins
        bins = bins.union(new_bins)

        # evaluate algorithm effectiveness
        bin_fraction_explored = len(bins) / num_bins ** 2
        if verbose:
            print_block('%s GENERATION %s: %5.2f%%' % (run_id, gen, bin_fraction_explored * 100))
        if bin_fraction_explored >= next_benchmark:
            print_block("%s: %5.2f%% exploration accomplished at generation %d" %
                ('{:%Y-%m-%d %H:%M:%S}'.format(datetime.now()), bin_fraction_explored * 100, gen))
            if benchmarks:
                next_benchmark = benchmarks.pop(0)
            else:
                last_benchmark_reached = True

        if config['output_all_graphs']:
            output_path = os.path.join(config['output_dir'], "binplot_%d.png" % gen)
            delaunay_figure(box_r, num_bins, output_path, children=new_box_r, parents=parents_r,
                            bins=bin_counts, new_bins=new_bins,
                            title="Generation %d: %d/%d (+%d) %5.2f%% (+%5.2f %%)" %
                                (gen, len(bins), num_bins ** 2, len(new_bins),
                                100*float(len(bins)) / num_bins ** 2, 100*float(len(new_bins)) / num_bins ** 2 ),
                            patches=None, prop1range=prop1range, prop2range=prop2range, \
                            perturbation_methods=perturbation_methods, showtriangulation=False)

        if gen <= 10 or (gen <=50 and gen % 10 == 0) or (gen <=500 and gen % 50 == 0) or \
            gen % 100 == 0 or last_benchmark_reached:
            output_path = os.path.join(config['output_dir'], "triplot_%d.png" % gen)
            delaunay_figure(box_r, num_bins, output_path, children=new_box_r, parents=parents_r,
                            bins=bin_counts, new_bins=new_bins,
                            title="Generation %d: %d/%d (+%d) %5.2f%% (+%5.2f %%)" %
                                (gen, len(bins), num_bins ** 2, len(new_bins),
                                100*float(len(bins)) / num_bins ** 2, 100*float(len(new_bins)) / num_bins ** 2 ),
                            patches=None, prop1range=prop1range, prop2range=prop2range, \
                            perturbation_methods=perturbation_methods)

        box_d = np.append(box_d, new_box_d, axis=0)
        box_r = np.append(box_r, new_box_r, axis=0)

        if last_benchmark_reached:
            break
