#!/usr/bin/env python3
import csv
import sys

import click
from sqlalchemy.orm import joinedload

from htsohm import db
from htsohm.db import Material
from htsohm.htsohm_serial import calc_bin

@click.command()
@click.argument('database-path', type=click.Path())
@click.option('--start-id', default=0, type=int)
def output_csv(database_path, start_id=0):
    db.init_database(db.get_sqlite_dbcs(database_path))
    session = db.get_session()
    mats = session.query(Material) \
        .options(joinedload("structure").joinedload("lennard_jones")) \
        .options(joinedload("gas_loading")) \
        .options(joinedload("void_fraction")) \
        .filter(Material.id >= start_id).all()

    f = csv.writer(sys.stdout, lineterminator="\n")
    f.writerow(["id", "parent_id",  "a", "b", "c", "sigma", "epsilon", "void_fraction", "void_fraction_geo", "absolute_volumetric_loading", "absolute_volumetric_loading_error"])
    for m in mats:
        f.writerow([m.id, m.parent_id, m.structure.a, m.structure.b, m.structure.c,
            m.structure.lennard_jones[0].sigma, m.structure.lennard_jones[0].epsilon,
            m.void_fraction[0].void_fraction, m.void_fraction[0].void_fraction_geo,
            m.gas_loading[0].absolute_volumetric_loading, m.gas_loading[0].absolute_volumetric_loading_error
        ])

@click.command()
@click.argument('csv-path', type=click.Path())
@click.option('-b', '--bin', nargs=4, type=click.Tuple([int, float, float, int]), multiple=True, help="column lower_bound upper_bound num_bins")
def csv_add_bin(csv_path, bin):
    with open(csv_path) as f:
        csv_in = csv.reader(f)
        csv_out = csv.writer(sys.stdout, lineterminator="\n")
        header = next(csv_in)

        bin_col_labels = ["bin%d" % col for col, _, _, _ in bin]
        csv_out.writerow(header + bin_col_labels + ["unique_bins"])
        unique_bins = set()
        for row in csv_in:
            calcd_bins = [calc_bin(float(row[col]), lb, ub, nb) for col, lb, ub, nb in bin]
            unique_bins = unique_bins.union(set([tuple(calcd_bins)]))
            csv_out.writerow(row + calcd_bins + [len(unique_bins)])


if __name__ == '__main__':
    output_csv()
