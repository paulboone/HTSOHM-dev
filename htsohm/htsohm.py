import sys

from htsohm import config
from htsohm.db import session, Material
from htsohm import pseudomaterial_generator
from htsohm.simulation.run_all import run_all_simulations

def evaluate_convergence(run_id):
    '''Determines convergence by calculating fraction of empty bins remaining.
    
    Args:
        run_id (str): identification string for run.

    Returns:
        bool: True if fraction of empty bins is less-than-or-equal-to the convergence cutoff
            specified in the configuration file.
    '''
    # query number of occupied bins for considering each of the simulation-types being run
    query_group = []
    for simulation in config["simulations"]:
        query_group.append(getattr(Material, '{}_bin'.format(simulation)))
    occupied_bins = session.query(*query_group).distinct() \
        .filter(Material.run_id == run_id).group_by(*query_group).count()

    # calculate fraction of empty bins remaining, compare to convergence cutoff
    total_bins = config['number_of_convergence_bins'] ** len(query_group)
    fraction_empty_bins = (total_bins - occupied_bins) / total_bins

    print("{0}\nFRACTION OF EMPTY BINS REMAINING : {1}\n{0}".format("=" * 80, fraction_empty_bins))

    return fraction_empty_bins <= config['convergence_cutoff_criteria']

def worker_run_loop(run_id):
    """
    Args:
        run_id (str): identification string for run.

    Manages overall routine for generating pseudomaterials and simulating their properties.
    Method runs until convergence cutoff is reached.

    """
    converged = False
    # run until convergence cutoff is reached
    while not converged:
        # generate pseudomaterial
        material, structure = pseudomaterial_generator.random.new_material(run_id,
                config["structure_parameters"])

        # simulate properties of interest
        run_all_simulations(material, structure)

        # add pseudomaterial to database
        session.add(material)
        session.commit()

        sys.stdout.flush()

        # evaluate convergence
        converged = evaluate_convergence(run_id)
