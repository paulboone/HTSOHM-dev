from glob import glob
import os
from uuid import uuid4

from htsohm import db, load_config_file
from htsohm.simulation import void_fraction


config_path = "../settings-research/ml_vf_best_alldof_small.yaml"
material_id = 23975

config = load_config_file(config_path)
db.init_database(db.get_sqlite_dbcs())
session = db.get_session()

from htsohm.db import Material

m = session.query(Material).get(material_id)


output_dir = "output_{}_{}".format(m.uuid, uuid4())
os.makedirs(output_dir, exist_ok=True)
void_fraction.write_output_files(m, config["simulations"][1], output_dir)
