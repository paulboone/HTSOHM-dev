#!/usr/bin/env python3

import csv
from itertools import chain
import sys

import click
import numpy as np

from htsohm import load_config_file, db
from htsohm.db import Material
from htsohm.htsohm_serial import calc_bins


def dof_analysis(config_path, run_id=None):
    config = load_config_file(config_path)
    db.init_database(config["database_connection_string"])
    session = db.get_session()

    children_per_generation = config['children_per_generation']
    prop1range = config['prop1range']
    prop2range = config['prop2range']

    num_bins = config['number_of_convergence_bins']
    bin_counts = np.zeros((num_bins, num_bins))

    materials = session.query(Material)
    if run_id:
        materials = materials.filter(Material.run_id==run_id)

    mats = materials.all()
    mats_d = mats[0:children_per_generation]
    # bins_d =

    perturbation_types = ["lattice", "atom_types", "atom_sites", "density", "all"]
    tsv = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
    tsv.writerow([""] + list(chain.from_iterable([[t] * 4 for t in perturbation_types])))
    tsv.writerow(["gen"] + list(chain.from_iterable([["#", "∆vf", "∆ml", "dist"] for t in perturbation_types])))

    gen = 1
    new_mats_d = mats[gen*children_per_generation:(gen + 1)*children_per_generation]
    while len(new_mats_d) > 0:
        # num_materials, ∆vf, ∆ml, ∆all
        gen_stats = {t:[0, 0.0, 0.0, 0.0] for t in perturbation_types}
        for m in new_mats_d:
            m_stats = gen_stats[m.perturbation]
            m_stats[0] += 1
            dvf = abs(m.void_fraction[0].void_fraction - m.parent.void_fraction[0].void_fraction)
            dml = abs(m.gas_loading[0].absolute_volumetric_loading - m.parent.gas_loading[0].absolute_volumetric_loading)
            m_stats[1] += dvf
            m_stats[2] += dml
            m_stats[3] += (dvf ** 2 + dml ** 2) ** 0.5

        row = [gen] + list(chain.from_iterable([gen_stats[t] for t in perturbation_types]))
        tsv.writerow(row)

        gen += 1
        new_mats_d = mats[gen*children_per_generation:(gen + 1)*children_per_generation]


@click.command()
@click.argument('config_path', type=click.Path())
@click.option('--run_id', type=click.STRING, default=None)
def dof(config_path, run_id):
    dof_analysis(config_path, run_id)

if __name__ == '__main__':
    dof()
