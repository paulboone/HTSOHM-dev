#!/usr/bin/env python3

import csv
from itertools import chain
import os
import sys

import click
import numpy as np

from htsohm import load_config_file, db
from htsohm.db import Material
from htsohm.htsohm_run import calc_bins

def dof_analysis(config_path, output_directory):
    config = load_config_file(config_path)
    db.init_database(config["database_connection_string"], config["properties"])
    session = db.get_session()

    children_per_generation = config['children_per_generation']
    prop1range = config['prop1range']
    prop2range = config['prop2range']

    num_bins = config['num_bins']
    bin_counts = np.zeros((num_bins, num_bins))

    vf_binunits = (prop1range[1] - prop1range[0]) / num_bins
    ml_binunits = (prop2range[1] - prop2range[0]) / num_bins

    materials = session.query(Material)

    perturbation_types = ["lattice", "lattice_nodens", "atom_types", "atom_sites", "density", "all"]

    tsv_output_path = os.path.join(output_directory, "data.tsv")
    tsvfile = open(tsv_output_path, 'w')
    tsv = csv.writer(tsvfile, delimiter="\t", lineterminator="\n")
    tsv.writerow([""] + list(chain.from_iterable([[t] * 5 for t in perturbation_types])))
    tsv.writerow(["gen"] + list(chain.from_iterable([["#", "∆vf", "∆ml", "dist", "new_bins"] for t in perturbation_types])))

    ids = materials.all()
    props = [(m.void_fraction[0].void_fraction, m.gas_loading[0].absolute_volumetric_loading) for m in ids]

    new_ids = ids[0:children_per_generation]
    new_props = props[0:children_per_generation]
    new_bins = calc_bins(new_props, num_bins, prop1range=prop1range, prop2range=prop2range)
    for i, (bx, by) in enumerate(new_bins):
        bin_counts[bx,by] += 1

    pts = {t:[] for t in perturbation_types}
    gen = 1
    new_ids = ids[gen*children_per_generation:(gen + 1)*children_per_generation]
    new_props = props[gen*children_per_generation:(gen + 1)*children_per_generation]
    animation = [[[b[0], b[1], -1, -1] for b in new_bins]]

    while len(new_ids) > 0:
        new_bins = calc_bins(new_props, num_bins, prop1range=prop1range, prop2range=prop2range)

        gen_animation = []
        gen_stats = {t:[0, 0.0, 0.0, 0.0, 0] for t in perturbation_types}
        for i, m in enumerate(new_ids):
            m_stats = gen_stats[m.perturbation]
            m_stats[0] += 1
            dvf = (m.void_fraction[0].void_fraction - m.parent.void_fraction[0].void_fraction) / vf_binunits
            dml = (m.gas_loading[0].absolute_volumetric_loading - m.parent.gas_loading[0].absolute_volumetric_loading) / ml_binunits
            m_stats[1] += dvf
            m_stats[2] += dml
            m_stats[3] += (dvf ** 2 + dml ** 2) ** 0.5
            if bin_counts[new_bins[i][0], new_bins[i][1]] == 0:
                m_stats[4] += 1

            # generate information for animation script
            parent_r = (m.parent.void_fraction[0].void_fraction, m.parent.gas_loading[0].absolute_volumetric_loading)
            parent_bin = calc_bins([parent_r], num_bins, prop1range=prop1range, prop2range=prop2range)[0]
            gen_animation.append([new_bins[i][0], new_bins[i][1], parent_bin[0], parent_bin[1]])

            # this and dml needed for output of numpy arrays # num_materials, ∆vf, ∆ml, ∆all, new_bins
            pts[m.perturbation].append([m.parent.gas_loading[0].absolute_volumetric_loading / ml_binunits, dml])

        for i, (bx, by) in enumerate(new_bins):
            bin_counts[bx,by] += 1

        row = [gen] + list(chain.from_iterable([gen_stats[t] for t in perturbation_types]))
        tsv.writerow(row)

        gen += 1
        new_ids = ids[gen*children_per_generation:(gen + 1)*children_per_generation]
        new_props = props[gen*children_per_generation:(gen + 1)*children_per_generation]
        animation.append(gen_animation)

    np.save(os.path.join(output_directory, "animation"), animation)

    for k in pts:
        np.save(os.path.join(output_directory, k), pts[k])

@click.command()
@click.argument('config_path', type=click.Path())
@click.argument('output_directory', type=click.Path())
def dof(config_path, output_directory):
    dof_analysis(config_path, output_directory)

if __name__ == '__main__':
    dof()
