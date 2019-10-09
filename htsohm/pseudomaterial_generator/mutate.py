from math import isclose
from random import choice, random, uniform

import numpy as np
from scipy.spatial import Delaunay, distance
from sqlalchemy import text
from sqlalchemy.orm.session import make_transient

from htsohm import db
from htsohm.db import Material, Structure, LennardJones, AtomSite
from htsohm.pseudomaterial_generator.utilities import random_number_density

def random_position(x0, x1, strength):
    # get minimum distance between two points of (1) within the box, and (2) across the box boundary
    dx = min(abs(x0 - x1), 1 - abs(x0 - x1))
    if x0 > x1 and (x0 - x1) > 0.5: # then dx will be through boundary, so move right
        x2 = (x0 + strength * dx) % 1.
    if x0 >= x1 and (x0 - x1) < 0.5: # then dx will be in box, so move left
        x2 = x0 - strength * dx
    if x0 < x1 and (x1 - x0) > 0.5: # then dx will be through boundary, so move left
        x2 = (x0 - strength * dx) % 1.
    if x0 < x1 and (x1 - x0) < 0.5: # then dx will be in box, so move right
        x2 = x0 + strength * dx
    return x2

def net_charge(atom_sites):
    return sum([e.q for e in atom_sites])

def perturb_unweighted(curr_val, max_change, var_limits):
    new_val = curr_val + uniform(-max_change, max_change)
    return min(max(new_val, var_limits[0]), var_limits[1])


def mutate_material(run_id, parent_id, config):
    """Create records for pseudomaterial simulation and structure data."

    Args:
        run_id (str): identification string for run.
        number_of_atomtypes (int): number of different chemical species used to
            populate the unit cell.

    Returns:
        material (sqlalchemy.orm.query.Query): database row for storing
            simulation data specific to the material. See
            `htsohm/db/material.py` for more information.

    """
    number_of_atom_types    = config["number_of_atom_types"]
    lattice_limits          = config["lattice_constant_limits"]
    number_density_limits   = config["number_density_limits"]
    epsilon_limits          = config["epsilon_limits"]
    sigma_limits            = config["sigma_limits"]
    max_charge              = config["charge_limit"]
    strength                = config["mutation_strength"]
    perturb                 = config["perturb"]

    ########################################################################

    # get parent structure
    session = db.get_session()
    parent = session.query(Material).get(int(parent_id))
    ps = parent.structure

    # create database row
    child = parent.clone()
    cs = child.structure

    perturb = set(config["perturb"])

    if config["perturb_type"] == "random":
        child.perturbation = choice(perturb)
        perturb = {child.perturbation}
    else:
        child.perturbation = "all"

    print("Parent ID: %d" % (child.parent_id))
    print("PERTURBING: %s [%s]" % (child.perturbation, perturb))

    # perturb lattice constants
    if perturb & {"lattice", "lattice_nodens"}:
        cs.a = perturb_unweighted(cs.a, strength * (lattice_limits[1] - lattice_limits[0]), lattice_limits)
        cs.b = perturb_unweighted(cs.b, strength * (lattice_limits[1] - lattice_limits[0]), lattice_limits)
        cs.c = perturb_unweighted(cs.c, strength * (lattice_limits[1] - lattice_limits[0]), lattice_limits)

        if "fix_atoms" not in config and perturb & {"lattice"}:
            new_density = len(cs.atom_sites) / cs.volume
            child.number_density = min(max(new_density, number_density_limits[0]), number_density_limits[1])


    # store unit cell volume to row
    child.unit_cell_volume = cs.volume

    # perturb lennard-jones parameters
    if perturb & {"atom_types"}:
        for at in cs.lennard_jones:
            at.sigma = perturb_unweighted(at.sigma, strength * (sigma_limits[1] - sigma_limits[0]), sigma_limits)
            at.epsilon = perturb_unweighted(at.epsilon, strength * (epsilon_limits[1] - epsilon_limits[0]), epsilon_limits)


    # adjust # of atom sites to match density--should only be required if number density is perturbed!
    if "fix_atoms" in config:
        number_of_atoms = config['fix_atoms']
    else:
        # perturb number density/ number of atom-sites
        if perturb & ["density"]:
            child.number_density = perturb_unweighted(parent.number_density, \
                                    (number_density_limits[1] - number_density_limits[0])*strength, \
                                    number_density_limits)

        number_of_atoms = max(1, round(child.number_density * child.unit_cell_volume))

    if number_of_atoms < len(cs.atom_sites):
        print("***** deleting atom sites")
        cs.atom_sites = np.random.choice(cs.atom_sites, number_of_atoms, replace=False).tolist()
    elif number_of_atoms > len(cs.atom_sites):
        print("***** adding atom sites")
        for i in range(number_of_atoms - len(cs.atom_sites)):
            # TODO: new points always have ZERO charge?
            atom_type = "A_{}".format(choice(range(number_of_atom_types)))
            cs.atom_sites.append(AtomSite(atom_type=atom_type, x=random(), y=random(), z=random(), q=0.,
                lennard_jones=cs.get_lennard_jones(atom_type)))

    # remove atom-sites, if necessary
    # adjust charges if atom-sites were removed
    # while not isclose(net_charge(cs.atom_sites), 0., abs_tol=1.0e-6):
    #         # pick an atom-site at random
    #         i = choice(range(len(cs.atom_sites)))
    #         if net_charge(cs.atom_sites) < 0:
    #             dq = float("{0:.6f}".format(min(max_charge - cs.atom_sites[i].q, -net_charge(cs.atom_sites))))
    #             cs.atom_sites[i].q += dq
    #         else:
    #             dq = float("{0:.6f}".format(min(max_charge + cs.atom_sites[i].q, net_charge(cs.atom_sites))))
    #             cs.atom_sites[i].q -= dq

    # perturb atom-site positions
    if perturb & {"atom_sites"}:
        for a in cs.atom_sites:
            a.x = random_position(a.x, random(), strength)
            a.y = random_position(a.y, random(), strength)
            a.z = random_position(a.z, random(), strength)


    # calculate avg. sigma/epsilon values
    sigma_sum, epsilon_sum = 0, 0
    for a in cs.atom_sites:
        index = int(a.atom_type[2:])
        sigma_sum += cs.lennard_jones[index].sigma
        epsilon_sum += cs.lennard_jones[index].epsilon
    child.average_sigma = sigma_sum / number_of_atoms
    child.average_epsilon = epsilon_sum / number_of_atoms

    # TODO: USE NON-WEIGHTED PERTURBATION
    # perturb partial charges
    # for i in range(number_of_atoms):
    #     while True:
    #         try:
    #             random_q = uniform(-max_charge, max_charge)
    #             dq = strength * (random_q - cs.atom_sites[i].q)
    #             j = choice(range(number_of_atoms))
    #             # TODO: what if i == j ?
    #             if abs(cs.atom_sites[i].q + dq) <= max_charge and abs(cs.atom_sites[j].q - dq) <= max_charge:
    #                 cs.atom_sites[i].q += dq
    #                 cs.atom_sites[j].q -= dq
    #                 break
    #         except:
    #             pass

    print("PARENT UUID :\t{}".format(parent.uuid))
    print("CHILD UUID  :\t{}".format(child.uuid))
    print("lattice constants: (%.2f, %.2f, %.2f) => (%.2f, %.2f, %.2f)" % (ps.a, ps.b, ps.c, cs.a, cs.b, cs.c))
    print("number_density: %.2e => %.2e" % (parent.number_density, child.number_density))
    print("number of atoms: %.2f => %.2f" % (int(parent.number_density * parent.unit_cell_volume), number_of_atoms))
    parent_ljs = ", ".join(["(%.1f, %.1f)" % (ljs.epsilon, ljs.sigma) for ljs in ps.lennard_jones])
    child_ljs = ", ".join(["(%.1f, %.1f)" % (ljs.epsilon, ljs.sigma) for ljs in cs.lennard_jones])
    print("lennard jones: %s => %s" % (parent_ljs, child_ljs))

    # print("FRAMEWORK NET CHARGE :\t{}".format(sum([e.q for e in cs.atom_sites])))
    return child
