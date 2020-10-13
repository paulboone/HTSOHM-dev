#!/usr/bin/env python3
from math import log10

import click
import numpy as np

from htsohm import load_config_file, db
from htsohm.db import Material
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

    prop1range = config['prop1range']
    prop2range = config['prop2range']
    num_bins = config['num_bins']

    # vf_binunits = (prop1range[1] - prop1range[0]) / num_bins
    # ml_binunits = (prop2range[1] - prop2range[0]) / num_bins

    print("loading materials...")

    if csv_path:
        props = np.loadtxt(csv_path, delimiter=',', skiprows=1, usecols=(12,13,5,6), max_rows=last_material)
        print("%d rows loaded from csv" % props.shape[0])
        if sigma_limits:
            props = props[(sigma_limits[0] <= props[:,2]) & (props[:,2] <= sigma_limits[1])]
            print("%d rows after applying sigma limits" % props.shape[0])
        if epsilon_limits:
            props = props[(epsilon_limits[0] <= props[:,3]) & (props[:,3] <= epsilon_limits[1])]
            print("%d rows after applying epsilon limits" % props.shape[0])
    else:
        db.init_database(db.get_sqlite_dbcs(database_path), config["properties"])
        session = db.get_session()

        ids = session.query(Material).options(joinedload("void_fraction"), joinedload("henrys_coefficient"))
        if last_material:
            ids = ids.limit(last_material).all()
        else:
            ids = ids.all()

        print("calculating material properties...")
        props = [(m.void_fraction[0].get_void_fraction(),
            log10(max(m.henrys_coefficient[0].co2_henrys, 10.0 ** prop2range[0]))) for m in ids]


    last_generation_start = len(props) - last_children

    print("calculating bins...")
    bin_counts = np.zeros((num_bins, num_bins))
    start_bins = calc_bins(props[0:last_generation_start], num_bins, prop1range=prop1range, prop2range=prop2range)
    for i, (bx, by) in enumerate(start_bins):
        bin_counts[bx,by] += 1
    bins_explored = np.count_nonzero(bin_counts)
    new_bins = calc_bins(props[last_generation_start:], num_bins, prop1range=prop1range, prop2range=prop2range)
    print(len(new_bins), len(start_bins), len(set(new_bins) - set(start_bins)))
    new_bins = set(new_bins) - set(start_bins)
    print("bins explored = %d" % bins_explored)

    children = []
    parents = []
    if last_children > 0:
        children = np.array(props[last_generation_start:])
        parent_ids = np.array([m.parent_id for m in ids[last_generation_start:]])
        parents = np.array([props[pid - 1] for pid in parent_ids])

    addl_data = None
    if addl_data_path:
        print("adding additional data from: %s" % addl_data_path)
        addl_data = np.loadtxt(addl_data_path, delimiter=",", skiprows=1, usecols=(1,2))

    print("outputting graph...")
    output_path = "binplot_%d_materials.png" % len(props)
    delaunay_figure(props, num_bins, output_path, bins=bin_counts, new_bins=new_bins,
                    title="%d Materials: %d/%d %5.2f%%" % (len(props), bins_explored,
                    num_bins ** 2, 100*float(bins_explored / num_bins ** 2)),
                    prop1range=prop1range, prop2range=prop2range, show_triangulation=False, show_hull=False,
                    addl_data_set=addl_data, children=children, parents=parents)

if __name__ == '__main__':
    bin_graph()
