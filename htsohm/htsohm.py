import sys

from sqlalchemy.sql import func

from htsohm import config
from htsohm.db import session, Material
from htsohm import pseudomaterial_generator
from htsohm.simulation.run_all import run_all_simulations

def materials_in_generation(run_id, generation):
    """Count number of materials in a generation.

    Args:
        run_id (str): identification string for run.
        generation (int): iteration in overall bin-mutate-simulate rountine.

    Returns:
        Number(int) of materials in a particular generation that are present in
        the database (the final step in bin-mutate-simulate routine).

    """
    return session.query(Material) \
            .filter(Material.run_id == run_id,
                    Material.generation == generation).count()

def last_generation(run_id):
    """Finds latest generation present in database.

    Args:
        run_id (str): identification string for run.

    Returns:
        Last generation(int) to be included in database.

    """
    return session.query(func.max(Material.generation)) \
            .filter(Material.run_id == run_id,)[0][0]

def evaluate_convergence(run_id, generation):
    '''Determines convergence by calculating fraction of empty bins remaining.
    
    Args:
        run_id (str): identification string for run.
        generation (int): iteration in overall routine.

    Returns:
        bool: True if fraction of empty bins is less-than-or-equal-to the convergence cutoff
            specified in the configuration file.
    '''
    # query number of occupied bins for considering each of the simulation-types being run
    query_group = []
    for simulation in config["simulations"]:
        query_group.append(getattr(Material, '{}_bin'.format(simulation)))
    occupied_bins = session.query(*query_group).distinct() \
        .filter(Material.run_id == run_id, Material.generation < generation,
                Material.generation_index < config['children_per_generation']) \
        .group_by(*query_group).count()

    # calculate fraction of empty bins remaining, compare to convergence cutoff
    total_bins = config['number_of_convergence_bins'] ** len(query_group)
    fraction_empty_bins = (total_bins - occupied_bins) / total_bins
    return fraction_empty_bins <= config['convergence_cutoff_criteria']

def print_block(string):
    print('{0}\n{1}\n{0}'.format('=' * 80, string))

def worker_run_loop(run_id):
    """
    Args:
        run_id (str): identification string for run.

    Manages overall routine for generating pseudomaterials and simulating their properties.
    Method runs until convergence cutoff or maximum number of generations is reached.

    """
    print('CONFIG\n{0}'.format(config))

    gen = last_generation(run_id) or 0

    converged = False
    while not converged:
        print_block('GENERATION {}'.format(gen))
        size_of_generation = config['children_per_generation']

        while materials_in_generation(run_id, gen) < size_of_generation:

            material, structure = pseudomaterial_generator.random.new_material(run_id, gen, config["structure_parameters"])

            # simulate material properties
            run_all_simulations(material, structure)
            
            session.add(material)
            session.commit()

            # add material to database, as needed
            material.generation_index = material.calculate_generation_index()
            if material.generation_index < config['children_per_generation']:
                print_block('ADDING MATERIAL {}'.format(material.uuid))
                session.add(material)

            else:
                # delete excess rows
                # session.delete(material)
                pass
            session.commit()
            sys.stdout.flush()
        gen += 1
        converged = evaluate_convergence(run_id, gen)
