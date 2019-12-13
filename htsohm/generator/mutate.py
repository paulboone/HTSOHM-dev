from random import choice, random, uniform

from htsohm.slog import slog
from htsohm.generator.random import random_atom_sites, random_atom_types

def random_position(x0, x1, mutation_strength):
    # get minimum distance between two points of (1) within the box, and (2) across the box boundary
    dx = min(abs(x0 - x1), 1 - abs(x0 - x1))
    if x0 > x1 and (x0 - x1) > 0.5: # then dx will be through boundary, so move right
        x2 = (x0 + mutation_strength * dx) % 1.
    if x0 >= x1 and (x0 - x1) <= 0.5: # then dx will be in box, so move left
        x2 = x0 - mutation_strength * dx
    if x0 < x1 and (x1 - x0) >= 0.5: # then dx will be through boundary, so move left
        x2 = (x0 - mutation_strength * dx) % 1.
    if x0 < x1 and (x1 - x0) < 0.5: # then dx will be in box, so move right
        x2 = x0 + mutation_strength * dx
    return x2

def net_charge(atom_sites):
    return sum([e.q for e in atom_sites])

def perturb_unweighted(curr_val, max_change, var_limits):
    new_val = curr_val + uniform(-max_change, max_change)
    return min(max(new_val, var_limits[0]), var_limits[1])

def print_parent_child_diff(parent, child):
    ps = parent.structure
    cs = child.structure
    slog("PARENT UUID :\t{}".format(parent.uuid))
    slog("CHILD UUID  :\t{}".format(child.uuid))
    slog("lattice constants: (%.2f, %.2f, %.2f) => (%.2f, %.2f, %.2f)" % (ps.a, ps.b, ps.c, cs.a, cs.b, cs.c))
    slog("number of atoms: %.2f => %.2f" % (len(parent.structure.atom_sites),
                                             len(child.structure.atom_sites)))
    parent_ats = ", ".join(["(%.1f, %.1f)" % (ats.epsilon, ats.sigma) for ats in ps.atom_types])
    child_ats = ", ".join(["(%.1f, %.1f)" % (ats.epsilon, ats.sigma) for ats in cs.atom_types])
    slog("atom types: %s => %s" % (parent_ats, child_ats))

def mutate_material(parent, config):
    child = parent.clone()
    cs = child.structure

    perturb = set(config["perturb"])
    if config["perturb_type"] == "random":
        child.perturbation = choice(perturb)
        perturb = {child.perturbation}
    else:
        child.perturbation = "all"

    slog("Parent id: %d" % (child.parent_id))
    slog("Perturbing: %s [%s]" % (child.perturbation, perturb))
    ms = config["mutation_strength"]

    if config["number_of_atom_types"] > len(cs.atom_types):
        num_atom_types_to_add = config["number_of_atom_types"] - len(cs.atom_types)
        slog("Adding %d random atom types so we have number defined in the config" % num_atom_types_to_add)
        cs.atom_types += random_atom_types(num_atom_types_to_add, config)

    if perturb & {"num_atoms"} and random() < ms:
        if random() < 0.5: # remove an atoms
            if len(cs.atom_sites) > config['num_atoms_limits'][0]:
                site_to_remove = choice(cs.atom_sites)
                slog("Removing atom site: ", site_to_remove)
                removed_site = cs.atom_sites.remove(site_to_remove)
        else: # add an atom
            if len(cs.atom_sites) < config['num_atoms_limits'][1]:
                slog("Adding atom site...")
                cs.atom_sites += random_atom_sites(1, cs.atom_types)


    if perturb & {"atom_type_assignments"}:
        for i, atom in enumerate(cs.atom_sites):
            if random() < ms**2:
                new_atom_type = choice(cs.atom_types)
                slog("Reassigning atom type for site %d from %d to %d" %
                    (i, cs.atom_types.index(atom.atom_types), cs.atom_types.index(new_atom_type)))
                atom.atom_types = new_atom_type

    if perturb & {"atom_types"}:
        sigl = config["sigma_limits"]
        epsl = config["epsilon_limits"]
        for at in cs.atom_types:
            at.sigma = perturb_unweighted(at.sigma, ms * (sigl[1] - sigl[0]), sigl)
            at.epsilon = perturb_unweighted(at.epsilon, ms * (epsl[1] - epsl[0]), epsl)

    if perturb & {"lattice"}:
        ll = config["lattice_constant_limits"]
        cs.a = perturb_unweighted(cs.a, ms * (ll[1] - ll[0]), ll)
        cs.b = perturb_unweighted(cs.b, ms * (ll[1] - ll[0]), ll)
        cs.c = perturb_unweighted(cs.c, ms * (ll[1] - ll[0]), ll)
        child.number_density = len(cs.atom_sites) / cs.volume

    if perturb & {"atom_sites"}:
        for a in cs.atom_sites:
            a.x = random_position(a.x, random(), ms)
            a.y = random_position(a.y, random(), ms)
            a.z = random_position(a.z, random(), ms)

    print_parent_child_diff(parent, child)
    return child
