#!/usr/bin/env python3
import csv
import sys

import click
from sqlalchemy.orm import joinedload

from htsohm import db
from htsohm.bins import calc_bin
from htsohm.db import Material, AtomSite

@click.command()
@click.argument('database-path', type=click.Path())
@click.option('--start-id', default=0, type=int)
def output_csv(database_path, start_id=0):
    db.init_database(db.get_sqlite_dbcs(database_path))
    session = db.get_session()
    output_csv_from_db(session, start_id)

def output_csv_from_db(session, start_id=0, output_file=sys.stdout):
    mats = session.query(Material) \
        .options(joinedload("structure").joinedload("atom_types", "atom_sites")) \
        .options(joinedload("gas_loading")) \
        .options(joinedload("void_fraction")) \
        .filter(Material.id >= start_id).all()

    f = csv.writer(output_file, lineterminator="\n")
    f.writerow(["id", "parent_id", "generation", "a", "b", "c", "volume", "atom_sites", "number_density",
                "total_epsilon", "epsilon_density", "void_fraction", "void_fraction_geo",
                "absolute_volumetric_loading", "absolute_volumetric_loading_error", "site_distribution"])
    for m in mats:
        f.writerow([m.id, m.parent_id, m.generation, m.structure.a, m.structure.b, m.structure.c, m.structure.volume,
            len(m.structure.atom_sites),  m.structure.number_density, m.structure.total_epsilon, m.structure.epsilon_density,
            m.void_fraction[0].void_fraction, m.void_fraction[0].void_fraction_geo,
            m.gas_loading[0].absolute_volumetric_loading, m.gas_loading[0].absolute_volumetric_loading_error,
            m.structure.site_distribution
        ])


@click.command()
@click.argument('database-path', type=click.Path())
def output_atom_sites_csv(database_path):
    db.init_database(db.get_sqlite_dbcs(database_path))
    session = db.get_session()
    output_atom_sites_csv_from_db(session)

def output_atom_sites_csv_from_db(session, output_file=sys.stdout):
    sites = session.query(AtomSite) \
        .options(joinedload("atom_types"))

    f = csv.writer(output_file, lineterminator="\n")
    f.writerow(["id", "structure_id", "x", "y", "z", "epsilon", "sigma", "a"])
    for s in sites:
        f.writerow([s.id, s.structure_id, s.x, s.y, s.z, s.atom_types.epsilon, s.atom_types.sigma, s.structure.a])

@click.command()
@click.argument('csv-path', type=click.Path())
@click.option('-b', '--bin', nargs=4, type=click.Tuple([int, float, float, int]), multiple=True, help="column lower_bound upper_bound num_bins")
def csv_add_bin(csv_path, bin, output_file=sys.stdout):
    csv_add_bin_column(csv_path, bin, output_file)

def csv_add_bin_column(csv_path, bin, output_file=sys.stdout):
    with open(csv_path) as f:
        csv_in = csv.reader(f)
        csv_out = csv.writer(output_file, lineterminator="\n")
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
