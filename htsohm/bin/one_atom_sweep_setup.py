#!/usr/bin/env python3

import click
import numpy as np
from sqlalchemy.orm import joinedload

from htsohm import load_config_file, db
from htsohm.db import Material

@click.command()
@click.argument('config-path', type=click.Path())
def sweep_materials(config_path):
    config = load_config_file(config_path)
    db.init_database(config["database_connection_string"])
    session = db.get_session()

    if 'sweep_points' in config:
        lattice_sweep_points = sigma_sweep_points = epsilon_sweep_points = config['sweep_points']
    else:
        lattice_sweep_points = config['lattice_sweep_points']
        sigma_sweep_points = config['sigma_sweep_points']
        epsilon_sweep_points = config['epsilon_sweep_points']

    scfg = config['structure_parameters']
    eps_d = np.linspace(*scfg['epsilon_limits'], epsilon_sweep_points)
    sig_d = np.linspace(*scfg['sigma_limits'], sigma_sweep_points)
    a_d = np.linspace(*scfg['lattice_constant_limits'], lattice_sweep_points)
    # always do symmetrical with one-atom only
    lattice_coords = [(a, a, a) for a in a_d]

    # b_d = np.linspace(*scfg['lattice_constant_limits'], config['sweep_points'])
    # c_d = np.linspace(*scfg['lattice_constant_limits'], config['sweep_points'])
    # lattice_coords = np.array(np.meshgrid(a_d, b_d, c_d)).T.reshape(-1,3)

    # # remove symmetrical points
    # lattice_coords = map(sorted, lattice_coords)
    # lattice_coords = set(map(tuple, lattice_coords))

    print("epsilons: ", eps_d)
    print("sigmas: ", sig_d)
    print("lattice constants: ", a_d)

    for eps in eps_d:
        for sig in sig_d:
            for coords in lattice_coords:
                material = Material.one_atom_new(sig, eps, *coords)
                session.add(material)

    session.commit()

if __name__ == '__main__':
    sweep_materials()
