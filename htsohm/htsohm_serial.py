
import os
from datetime import datetime
import math

import numpy as np

from htsohm import pseudomaterial_generator, load_config_file
from htsohm.db import session, Material
from htsohm.simulation.run_all import run_all_simulations

from htsohm.figures import delaunay_figure
from htsohm.select.triangulation import choose_parents


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

def serial_runloop(config_path):
    """
    Args:
        run_id (str): identification string for run.

    """

    config = load_config_file(config_path)
    num_bins = config['number_of_convergence_bins']
    run_id = datetime.now().isoformat()
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
    box_d = np.zeros(children_per_generation)
    box_r = -1 * np.ones((children_per_generation, 2))
    for i in range(children_per_generation):
        material = pseudomaterial_generator.random.new_material(run_id, config["structure_parameters"])
        run_all_simulations(material, config)
        session.add(material)
        session.commit()
        box_d[i] = material.id
        box_r[i,:] = (material.void_fraction[0].void_fraction, material.gas_loading[0].absolute_volumetric_loading)
        # box_r[i,:] = (material[prop1], material[prop2])

    bins = set(calc_bins(box_r, num_bins, prop1range=prop1range, prop2range=prop2range))
    print("bins", bins)

    last_benchmark_reached = False
    os.makedirs(config['visualization_output_dir'], exist_ok=True)

    for gen in range(1, config['max_generations'] + 1):

        parents_d, parents_r = choose_parents(children_per_generation, box_d, box_r, config['simplices_or_hull'])

        # mutate materials and simulate properties
        new_box_d = np.zeros(children_per_generation)
        new_box_r = -1 * np.ones((children_per_generation, 2))
        perturbation_methods = [""] * children_per_generation
        for i in range(children_per_generation):
            material = pseudomaterial_generator.mutate.mutate_material(run_id, parents_d[i], config["structure_parameters"])
            perturbation_methods[i] = material.perturbation
            run_all_simulations(material, config)
            session.add(material)
            session.commit()
            new_box_d[i] = material.id
            new_box_r[i,:] = (material.void_fraction[0].void_fraction, material.gas_loading[0].absolute_volumetric_loading)
            # new_box_r[i,:] = (material[prop1], material[prop2])

        # TODO: bins for methane loading?
        new_bins = set(calc_bins(new_box_r, num_bins, prop1range=prop1range, prop2range=prop2range)) - bins
        bins = bins.union(new_bins)

        # evaluate algorithm effectiveness
        bin_count = len(bins)
        bin_fraction_explored = bin_count / num_bins ** 2
        if verbose:
            print_block('%s GENERATION %s: %5.2f%%' % (run_id, gen, bin_fraction_explored * 100))
        if bin_fraction_explored >= next_benchmark:
            print_block("%s: %5.2f%% exploration accomplished at generation %d" %
                ('{:%Y-%m-%d %H:%M:%S}'.format(datetime.now()), bin_fraction_explored * 100, gen))
            if benchmarks:
                next_benchmark = benchmarks.pop(0)
            else:
                last_benchmark_reached = True

        if gen <= 10 or (gen <=50 and gen % 10 == 0) or (gen <=500 and gen % 50 == 0) or \
            gen % 100 == 0 or last_benchmark_reached:
            output_path = os.path.join(config['visualization_output_dir'], "triplot_%d.png" % gen)
            delaunay_figure(box_r, num_bins, output_path, children=new_box_r, parents=parents_r,
                            bins=bins, new_bins=new_bins,
                            title="Generation %d: %d/%d (+%d) %5.2f%% (+%5.2f %%)" %
                                (gen, len(bins), num_bins ** 2, len(new_bins),
                                100*float(len(bins)) / num_bins ** 2, 100*float(len(new_bins)) / num_bins ** 2 ),
                            patches=None, prop1range=prop1range, prop2range=prop2range, \
                            perturbation_methods=perturbation_methods)

        box_d = np.append(box_d, new_box_d, axis=0)
        box_r = np.append(box_r, new_box_r, axis=0)

        if last_benchmark_reached:
            break
