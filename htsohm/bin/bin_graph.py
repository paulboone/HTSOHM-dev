#!/usr/bin/env python3

import click
import numpy as np

from htsohm import load_config_file, db
from htsohm.db import Material, VoidFraction
from htsohm.figures import delaunay_figure
from htsohm.htsohm_run import calc_bins

from sqlalchemy.orm import joinedload

@click.command()
@click.argument('config-path', type=click.Path())
@click.option('--database-path', type=click.Path())
@click.option('--csv-path', type=click.Path())
@click.option('--last-material', '-l', type=int, default=None)
@click.option('--sigma-limits', type=float, nargs=2)
@click.option('--epsilon-limits', type=float, nargs=2)
@click.option('--addl-data-path', type=click.Path())
@click.option('--last-children', type=int, default=0)
def bin_graph(config_path, database_path=None, csv_path=None, last_material=None,
        sigma_limits=None, epsilon_limits=None, addl_data_path=None, last_children=0):

    config = load_config_file(config_path)
    VoidFraction.set_column_for_void_fraction(config['void_fraction_subtype'])

    prop1range = config['prop1range']
    prop2range = config['prop2range']
    num_bins = config['number_of_convergence_bins']

    # vf_binunits = (prop1range[1] - prop1range[0]) / num_bins
    # ml_binunits = (prop2range[1] - prop2range[0]) / num_bins

    print("loading materials...")

    if csv_path:
        mats_r = np.loadtxt(csv_path, delimiter=',', skiprows=1, usecols=(12,13,5,6), max_rows=last_material)
        print("%d rows loaded from csv" % mats_r.shape[0])
        if sigma_limits:
            mats_r = mats_r[(sigma_limits[0] <= mats_r[:,2]) & (mats_r[:,2] <= sigma_limits[1])]
            print("%d rows after applying sigma limits" % mats_r.shape[0])
        if epsilon_limits:
            mats_r = mats_r[(epsilon_limits[0] <= mats_r[:,3]) & (mats_r[:,3] <= epsilon_limits[1])]
            print("%d rows after applying epsilon limits" % mats_r.shape[0])
    else:
        db.init_database(db.get_sqlite_dbcs(database_path))
        session = db.get_session()

        mats_d = session.query(Material).options(joinedload("void_fraction"), joinedload("gas_loading"))
        if last_material:
            mats_d = mats_d.limit(last_material).all()
        else:
            mats_d = mats_d.all()

        print("calculating material properties...")
        mats_r = [(m.void_fraction[0].get_void_fraction(), m.gas_loading[0].absolute_volumetric_loading) for m in mats_d]


    last_generation_start = len(mats_r) - last_children

    print("calculating bins...")
    bin_counts = np.zeros((num_bins, num_bins))
    start_bins = calc_bins(mats_r[0:last_generation_start], num_bins, prop1range=prop1range, prop2range=prop2range)
    for i, (bx, by) in enumerate(start_bins):
        bin_counts[bx,by] += 1
    bins_explored = np.count_nonzero(bin_counts)
    new_bins = calc_bins(mats_r[last_generation_start:], num_bins, prop1range=prop1range, prop2range=prop2range)
    print(len(new_bins), len(start_bins), len(set(new_bins) - set(start_bins)))
    new_bins = set(new_bins) - set(start_bins)
    print("bins explored = %d" % bins_explored)

    children = []
    parents = []
    if last_children > 0:
        children = np.array(mats_r[last_generation_start:])
        parent_ids = np.array([m.parent_id for m in mats_d[last_generation_start:]])
        parents = np.array([mats_r[pid - 1] for pid in parent_ids])

    addl_data = None
    if addl_data_path:
        print("adding additional data from: %s" % addl_data_path)
        addl_data = np.loadtxt(addl_data_path, delimiter=",", skiprows=1, usecols=(1,2))

    print("outputting graph...")
    output_path = "binplot_%d_materials.png" % len(mats_r)
    delaunay_figure(mats_r, num_bins, output_path, bins=bin_counts, new_bins=new_bins,
                    title="%d Materials: %d/%d %5.2f%%" % (len(mats_r), bins_explored,
                    num_bins ** 2, 100*float(bins_explored / num_bins ** 2)),
                    prop1range=prop1range, prop2range=prop2range, show_triangulation=False, show_hull=False,
                    addl_data_set=addl_data, children=children, parents=parents)

if __name__ == '__main__':
    bin_graph()
