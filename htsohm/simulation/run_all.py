
from datetime import datetime

from htsohm.simulation import simulate
from htsohm.slog import slog

def run_all_simulations(material, config):
    """Simulate helium void fraction, gas loading, and surface area.

    Args:
        material (sqlalchemy.orm.query.Query): material to be analyzed.

    Depending on properties specified in config, adds simulated data for helium
    void fraction, gas loading, heat of adsorption, surface area, and
    corresponding bins to row in database corresponding to the input-material.
    """
    slog("-----------------------------------------------")
    for simulation_config in config["properties"]:
        slog('Time             : {:%Y-%m-%d %H:%M:%S}'.format(datetime.now()))
        slog("Simulation type  : {}".format(simulation_config["type"]))
        getattr(simulate, simulation_config["type"]).run(material, simulation_config, config)
        slog("--")

    slog('{:%Y-%m-%d %H:%M:%S}'.format(datetime.now()))
