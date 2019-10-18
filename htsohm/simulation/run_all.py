
from datetime import datetime

from htsohm.simulation import simulate

def run_all_simulations(material, config):
    """Simulate helium void fraction, gas loading, and surface area.

    Args:
        material (sqlalchemy.orm.query.Query): material to be analyzed.

    Depending on properties specified in config, adds simulated data for helium
    void fraction, gas loading, heat of adsorption, surface area, and
    corresponding bins to row in database corresponding to the input-material.

    """
    for simulation_number in config["simulations"]:
        print('{:%Y-%m-%d %H:%M:%S}'.format(datetime.now()))
        simulation_config = config["simulations"][simulation_number]
        getattr(simulate, simulation_config["type"]).run(material, simulation_config, config)
        #material.update_from_dict(results)
    print('{:%Y-%m-%d %H:%M:%S}'.format(datetime.now()))
