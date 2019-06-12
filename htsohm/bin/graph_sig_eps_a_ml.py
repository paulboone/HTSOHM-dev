#!/usr/bin/env python3

import click
import matplotlib.pyplot as plt
import numpy as np
from sqlalchemy.orm import joinedload

from htsohm import load_config_file, db
from htsohm.db import Material

def limit_index(v, limits, sweep_points):
    return int(round((sweep_points - 1) * (v - limits[0]) / (limits[1] - limits[0])))

@click.command()
@click.argument('config-path', type=click.Path())
@click.option('--database-path', type=click.Path())
@click.option('--addl-data-path', type=click.Path())
def bin_graph(config_path, database_path=None, addl_data_path=None):
    config = load_config_file(config_path)
    db.init_database(db.get_sqlite_dbcs(database_path))
    session = db.get_session()

    num_bins = config['number_of_convergence_bins']
    prop1range = [2.5, 5.0]
    prop2range = config['prop2range']

    print("loading materials...")
    mats = session.query(Material) \
        .options(joinedload("structure").joinedload("lennard_jones")) \
        .options(joinedload("gas_loading")).all()

    print("calculating material properties...")
    mats_by_lj = {}
    for m in mats:
        lj = (m.structure.lennard_jones[0].sigma, m.structure.lennard_jones[0].epsilon)
        if lj not in mats_by_lj:
            mats_by_lj[lj] = []

        mats_by_lj[lj].append([m.structure.a, m.gas_loading[0].absolute_volumetric_loading])

    print("plotting...")
    fig = plt.figure(figsize=(12,12), tight_layout=True)
    ax = fig.add_subplot(1, 1, 1)
    ax.set_xlim(prop1range[0], prop1range[1])
    ax.set_ylim(prop2range[0], prop2range[1])
    ax.set_xlabel("Lattice constant [Å]")
    ax.set_ylabel("Methane Loading (V [STP]/V)")
    ax.set_yticks(prop2range[1] * np.array([0.0, 0.25, 0.5, 0.75, 1.0]))
    ax.set_yticks(prop2range[1] * np.array(range(0,num_bins + 1))/num_bins, minor=True)

    ax.set_xticks(prop1range[0] + (prop1range[1] - prop1range[0]) * np.array([0.0, 0.25, 0.5, 0.75, 1.0]))
    ax.set_xticks(prop1range[0] + (prop1range[1] - prop1range[0]) * np.array(range(0,num_bins + 1))/num_bins, minor=True)

    # if show_grid:
    ax.grid(linestyle='-', color='0.8', zorder=0)


    absolute_limits_a = np.linspace(prop1range[0], prop1range[1],100)
    ml_atoms_a3_stp = 2.69e-5
    absolute_limits_ml = [ (1/a**3) / ml_atoms_a3_stp for a in absolute_limits_a]
    ax.plot(absolute_limits_a, absolute_limits_ml, lw=3, linestyle="--", color="black", zorder=15, label="all sites filled")
    # print(mats_by_lj.keys())

    for (sig, eps), a_ml in mats_by_lj.items():
        a_ml = np.array(a_ml)
        sig_index = limit_index(sig, config['structure_parameters']['sigma_limits'], config['sweep_points'])
        eps_index = limit_index(eps, config['structure_parameters']['epsilon_limits'], config['sweep_points'])
        # print(sig, eps, marker_index)


        alpha = (eps_index + 1) / config['sweep_points']
        if eps_index + 1 == config['sweep_points']:
            label = "sigma = %3.2f" % sig
        else:
            label = None
        #  marker="$%d$" % sig_index, markersize=10,
        ax.plot(a_ml[:,0], a_ml[:,1], lw=3, zorder=20, alpha=alpha, label=label)

    ax.legend()
    ax.set_title("Methane loading vs lattice constant. Lines colored by sigma. \nLine transparency shows epsilon (no transparency is highest epsilon value; highest transparency is lowest epsilon value)")

    fig.savefig("sig_eps_a_ml.png")
    plt.close(fig)


if __name__ == '__main__':
    bin_graph()
