from glob import glob
import os

import click

from htsohm import db, load_config_file
from htsohm.simulation import simulate


@click.command()
@click.argument('config-path', type=click.Path())
@click.argument('material_ids', type=int, nargs=-1)
@click.option('--database-path', type=click.Path())
def output_config_files(config_path, material_ids, database_path=None):
    config = load_config_file(config_path)
    db.init_database(db.get_sqlite_dbcs(database_path), config["properties"])
    session = db.get_session()

    from htsohm.db import Material

    for m_id in material_ids:
        m = session.query(Material).get(m_id)

        for i in config["simulations"]:
            simcfg = config["simulations"][i]
            output_dir = "output_%d_%s_%s_%d" % (m.id, simcfg["type"], i)
            os.makedirs(output_dir, exist_ok=True)

            sim = getattr(simulate, simcfg["type"])
            sim.write_input_files(m, simcfg, output_dir)


if __name__ == '__main__':
    output_config_files()
