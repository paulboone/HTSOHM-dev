import sys

from sqlalchemy.sql import text

from htsohm import config
from htsohm.db import engine, session, Material
from htsohm import pseudomaterial_generator
from htsohm.simulation.run_all import run_all_simulations

def add_columns(query_string, count, flag):
    for i in range(count):
        query_string += """
        , {}{}.bin_value""".format(flag, i)
    return query_string

def get_table(flag):
    if flag == "g":
        return "gas_loadings"
    elif flag == "s":
        return "surface_areas"
    elif flag == "v":
        return "void_fractions"
    else:
        print("WARNING : UNEXPECTED OUTPUT TABLE FLAG")

def add_joins(query_string, count, flag):
    for i in range(count):
        query_string += """
        join {0} {1}{2} on m.id={1}{2}.material_id""".format(get_table(flag),
                flag, i)
    return query_string

def add_filters(query_string, simulation_config, simulation_type, flag):
    count = 0
    for simulation_id in simulation_config:
        simulation_details = simulation_config[simulation_id]
        if simulation_type == simulation_details["type"]:
            filter_line = """
            and {}{}.adsorbate='{}'""".format(flag, count,
                    simulation_details["adsorbate"])
            if "pressure" in simulation_details:
                filter_line += """ and {}{}.pressure={}""".format(flag, count,
                        simulation_details["pressure"])
            if "temperature" in simulation_details:
                filter_line += """ and {}{}.temperature={}""".format(flag, count,
                        simulation_details["temperature"])
            query_string += filter_line
            count += 1
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
    bin_dimension = gas_loading_count + surface_area_count + void_fraction_count

    query_string = """
    select count(bins) from
      (select distinct m.run_id"""

    # add columns for each bin-dimension
    query_string = add_columns(query_string, gas_loading_count, "g")
    query_string = add_columns(query_string, surface_area_count, "s")
    query_string = add_columns(query_string, void_fraction_count, "v")

    query_string += """
       from materials m"""

    query_string = add_joins(query_string, gas_loading_count, "g")
    query_string = add_joins(query_string, surface_area_count, "s")
    query_string = add_joins(query_string, void_fraction_count, "v")

    # filter by run_id
    query_string += """
       where m.run_id=:run_id"""

    if gas_loading_count > 0:
        query_string = add_filters(query_string, config["simulations"], "gas_loading", "g")
    if surface_area_count > 0:
        query_string = add_filters(query_string, config["simulations"], "surface_area", "s")
    if void_fraction_count > 0:
        query_string = add_filters(query_string, config["simulations"], "void_fraction", "v")

    query_string += """
       ) as bins;"""

    print(query_string)

    sql = text(query_string) 
    rows = engine.connect().execute(
        sql,
        run_id = run_id).fetchall()
    occupied_bins = rows[0][0]

    # calculate fraction of empty bins remaining, compare to convergence cutoff
    total_bins = config['number_of_convergence_bins'] ** bin_dimension
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
