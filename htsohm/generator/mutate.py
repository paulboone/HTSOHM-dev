from math import ceil, floor
from random import choice, choices, random, uniform

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
    # the / 2 is to make sure that if the max_change is the entire variable range (i.e. 100%
    # mutation strength), then if we are in the center of the range, this equivalent to generating
    # a new random number in the range.
    new_val = curr_val + uniform(-max_change, max_change) / 2
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
    slog("Starting site density: %7.6f" % (cs.number_density))
    ms = config["mutation_strength"]

    if perturb & {"lattice"}:
        ll = config["lattice_constant_limits"]
        if config["lattice_cubic"]:
            new_a = perturb_unweighted(cs.a, ms * (ll[1] - ll[0]), ll)
            slog("a %3.2f => %3.2f" % (cs.a, new_a))
            # new_a needs to be bounded s.t. the # atoms doesn't violate the density limits
            new_a = max(new_a, (len(cs.atom_sites) / config["density_limits"][1])**(1/3))
            new_a = min(new_a, (len(cs.atom_sites) / config["density_limits"][0])**(1/3))
            slog("a %3.2f => %3.2f (after enforcing density_limits)" % (cs.a, new_a))
            cs.a = new_a
            cs.b = new_a
            cs.c = new_a
        else:
            raise(Exception("change needs to be made to support non-cubic lattices..."))
            cs.b = perturb_unweighted(cs.b, ms * (ll[1] - ll[0]), ll)
            cs.c = perturb_unweighted(cs.c, ms * (ll[1] - ll[0]), ll)

    if perturb & {"num_atoms"}:
        if random() < ms:
            if random() < 0.5: # remove an atoms
                if len(cs.atom_sites) > config['num_atoms_limits'][0]:
                    site_to_remove = choice(cs.atom_sites)
                    slog("Removing atom site: ", site_to_remove)
                    cs.atom_sites.remove(site_to_remove)
            else: # add an atom
                if len(cs.atom_sites) < config['num_atoms_limits'][1]:
                    slog("Adding atom site...")
                    cs.atom_sites += random_atom_sites(1, cs.atom_types)

    elif perturb & {"density"}:
        if random() < config["density_frequency"]:
            slog("Attempting density change...")
            dl = config["density_limits"]
            number_density = perturb_unweighted(cs.number_density, ms * (dl[1] - dl[0]), dl)
            slog("site-density %7.6f => %7.6f" % (cs.number_density, number_density))
            number_of_atoms = round(number_density * cs.volume)
            slog("# atoms %d => %d" % (len(cs.atom_sites), number_of_atoms))

            # num atoms needs to be bounded s.t. any integer # is within the density limit range, not just outside it
            number_of_atoms = max(number_of_atoms, ceil(config["density_limits"][0] * cs.volume))
            number_of_atoms = min(number_of_atoms, floor(config["density_limits"][1] * cs.volume))

            slog("# atoms %d => %d (after enforcing density_limits on rounding to int)" % (len(cs.atom_sites), number_of_atoms))
            delta_atoms = number_of_atoms - len(cs.atom_sites)

            slog("delta_atoms: %d" % delta_atoms)
            if 'max_delta_atoms' in config:
                # additionally bound delta
                delta_atoms = max(delta_atoms, -config["max_delta_atoms"])
                delta_atoms = min(delta_atoms, config["max_delta_atoms"])
                slog("delta_atoms: %d (after max_delta_atoms)" % delta_atoms)

            if delta_atoms < 0: # remove atoms
                slog("Removing %d atom site(s)..." % delta_atoms)
                for site_to_remove in choices(cs.atom_sites, k=delta_atoms):
                    cs.atom_sites.remove(site_to_remove)
            elif delta_atoms > 0: # add atoms
                slog("Adding %d atom site(s)..." % delta_atoms)
                cs.atom_sites += random_atom_sites(delta_atoms, cs.atom_types)


    if config["number_of_atom_types"] > len(cs.atom_types):
        num_atom_types_to_add = config["number_of_atom_types"] - len(cs.atom_types)
        slog("Adding %d random atom types so we have number defined in the config" % num_atom_types_to_add)
        cs.atom_types += random_atom_types(num_atom_types_to_add, config)

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

    if perturb & {"atom_sites"}:
        for a in cs.atom_sites:
            a.x = random_position(a.x, random(), ms)
            a.y = random_position(a.y, random(), ms)
            a.z = random_position(a.z, random(), ms)

    print_parent_child_diff(parent, child)
    return child
