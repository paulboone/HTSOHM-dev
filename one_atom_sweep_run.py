#!/usr/bin/env python3

import click

from htsohm import load_config_file, db
from htsohm.db import Material
from htsohm.simulation.run_all import run_all_simulations

@click.command()
@click.argument('config-path', type=click.Path())
@click.argument('range-start', type=int)
@click.argument('range-stop', type=int)
def run_materials(config_path, range_start, range_stop):
    config = load_config_file(config_path)
    db.init_database(config["database_connection_string"])
    session = db.get_session()

    mats = session.query(Material).filter(Material.id >= range_start, Material.id < range_stop).all()

    print(len(mats))
    for m in mats:
        run_all_simulations(m, config)
        session.add(m)
        session.commit()

if __name__ == '__main__':
    run_materials()
