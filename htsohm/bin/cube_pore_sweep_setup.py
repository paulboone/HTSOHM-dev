#!/usr/bin/env python3

import click
import numpy as np
from sqlalchemy.orm import joinedload

from htsohm import load_config_file, db
from htsohm.db import Material

@click.command()
@click.argument('config-path', type=click.Path())
def sweep_setup(config_path):
    config = load_config_file(config_path)
    db.init_database(config["database_connection_string"])
    session = db.get_session()

    directional_atom_sweep_points = config['directional_atom_sweep_points']
    sigma_sweep_points = config['sigma_sweep_points']
    epsilon_sweep_points = config['epsilon_sweep_points']

    scfg = config['structure_parameters']
    eps_d = np.linspace(*scfg['epsilon_limits'], epsilon_sweep_points)
    sig_d = np.linspace(*scfg['sigma_limits'], sigma_sweep_points)
    atoms_d = np.linspace(*scfg['directional_atom_limits'], directional_atom_sweep_points, dtype=int)
    atom_diameter = scfg['atom_diameter']

    print("epsilons: ", eps_d)
    print("sigmas: ", sig_d)
    print("num_atoms: ", atoms_d)

    for eps in eps_d:
        for sig in sig_d:
            for num_atoms in atoms_d:
                material = Material.cube_pore_new(sig, eps, num_atoms, atom_diameter)
                session.add(material)

    session.commit()

if __name__ == '__main__':
    sweep_setup()
