import sys

from sqlalchemy import func

from htsohm import config
from htsohm.db import session, Material
from htsohm import pseudomaterial_generator
from htsohm.simulation.run_all import run_all_simulations

def count_number_of_materials(run_id):
    return session.query(func.count(Material.id)).filter(Material.run_id==run_id)[0][0]

def worker_run_loop(run_id):
    """
    Args:
        run_id (str): identification string for run.

    Manages overall routine for generating pseudomaterials and simulating their properties.
    Method runs until convergence cutoff is reached.

    """
    while count_number_of_materials(run_id) < 1000000:
        if count_number_of_materials(run_id) < config["seed_count"]:
            # generate pseudomaterial
            material = pseudomaterial_generator.random.new_material(run_id, config["structure_parameters"])
        else:
            # mutate pseudomaterial
            material = pseudomaterial_generator.mutate.new_material(run_id, config)

        # simulate properties of interest
        run_all_simulations(material)

        # add pseudomaterial to database
        session.add(material)
        session.commit()

        sys.stdout.flush()
