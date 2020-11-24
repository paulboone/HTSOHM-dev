#!/usr/bin/env python3

import click
import numpy as np

from htsohm import load_config_file, db
from htsohm.db import Material, VoidFraction
from htsohm.htsohm_run import load_restart_db

from sqlalchemy.orm import joinedload

@click.command()
@click.argument('config-path', type=click.Path())
@click.argument('database-path', type=click.Path())
@click.option('--generation', '-g', type=int, default=None)
@click.option('-o', '--output-path', type=click.Path(), default="-")
def num_materials_per_bin(config_path, database_path, generation=None, output_path="-"):
    """outputs materials per bin for cover visualization script in blender"""

    config = load_config_file(config_path)
    prop1range = config['prop1range']
    prop2range = config['prop2range']
    num_bins = config['number_of_convergence_bins']
    VoidFraction.set_column_for_void_fraction(config['void_fraction_subtype'])

    engine, session = db.init_database(config["database_connection_string"])

    generation = 500

    _, _, bin_counts, _, _, _ = load_restart_db(generation, num_bins, prop1range, prop2range, session)

    with click.open_file(output_path, 'w') as f:
        np.savetxt(f, bin_counts, "%d", delimiter=",")


if __name__ == '__main__':
    num_materials_per_bin()
