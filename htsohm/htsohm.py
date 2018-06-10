import sys

from sqlalchemy.sql import text

from htsohm import config
from htsohm.db import engine, session, Material
from htsohm import pseudomaterial_generator
from htsohm.simulation.run_all import run_all_simulations

def add_subquery(query_string, number_of_queries, subquery_string):
    for i in range(number_of_queries):
        query_string += subquery_string.format(i)
    return query_string

def evaluate_convergence(run_id):
    '''Determines convergence by calculating fraction of empty bins remaining.
    
    Args:
        run_id (str): identification string for run.

    Returns:
        bool: True if fraction of empty bins is less-than-or-equal-to the convergence cutoff
            specified in the configuration file.
    '''
    gas_loading_count, surface_area_count, void_fraction_count = 0, 0, 0
    for simulation_id in config["simulations"]:
        simulation_type = config["simulations"][simulation_id]["type"]
        if simulation_type == "gas_loading":
            gas_loading_count += 1
        elif simulation_type == "surface_area":
            surface_area_count += 1
        elif simulation_type == "void_fraction":
            void_fraction_count += 1
        else:
            print("WARNING : UNEXPECTED SIMULATION TYPE.")
    
    query_string = """
    select m.id
    """
    
    gas_loading_subquery = """
       , (select g{0}.bin_value from gas_loadings g{0}
           where g{0}.material_id = m.id
           order by g{0}.adsorbate, g{0}.temperature, g{0}.pressure
           limit 1 offset {0}
         ) as g_bin_value_{0}
    """       
    if gas_loading_count > 0:
        query_string = add_subquery(query_string, gas_loading_count,
                                   gas_loading_subquery)
    
    surface_area_subquery = """
       , (select s{0}.bin_value from surface_areas s{0}
           where s{0}.material_id = m.id
           order by s{0}.adsorbate
           limit 1 offset {0}
         ) as s_bin_value_{0}
    """
    if surface_area_count > 0:
        query_string = add_subquery(query_string, surface_area_count,
                                   surface_area_subquery)
    
    void_fraction_subquery = """
       , (select v{0}.bin_value from void_fractions v{0}
           where v{0}.material_id = m.id
           order by v{0}.adsorbate, v{0}.temperature
           limit 1 offset {0}
         ) as v_bin_value_{0}
    """
    if void_fraction_count > 0:
        query_string = add_subquery(query_string, void_fraction_count,
                                   void_fraction_subquery)
    
    query_string += """
      from materials m
    order by m.id;
    """
    
    print(query_string)
    
    sql = text(query_string) 
    rows = engine.connect().execute(
        sql,
        run_id = run_id).fetchall()
    bins = []
    for row in rows:
        bin_ = [row[1:]]
        if bin_ not in bins:
            bins.append(bin_)
    occupied_bins = len(bins)

    # calculate fraction of empty bins remaining, compare to convergence cutoff
    total_bins = config['number_of_convergence_bins'] ** (gas_loading_count + surface_area_count + void_fraction_count)
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
