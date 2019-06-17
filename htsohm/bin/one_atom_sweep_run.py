#!/usr/bin/env python3

import click

from htsohm import load_config_file, db
from htsohm.db import Material
from htsohm.simulation.run_all import run_all_simulations

@click.command()
@click.argument('config-path', type=click.Path())
@click.option('--workers', type=int, nargs=2, default=(1,1), help='NUM_WORKERS THIS_WORKER_INDEX')
def run_materials(config_path, workers=(1,1)):
    config = load_config_file(config_path)
    db.init_database(config["database_connection_string"])
    session = db.get_session()

    mats = session.query(Material).all()
    num_workers, worker_num = workers
    mats =  mats[worker_num - 1::num_workers]

    print(len(mats))
    for m in mats:
        print("---------------")
        print("%d" % m.id)
        run_all_simulations(m, config)
        session.add(m)
        session.commit()

if __name__ == '__main__':
    run_materials()
