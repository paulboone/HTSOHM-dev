import math
from random import choice, random, uniform

from htsohm.slog import slog
from htsohm.generator.random import random_atom_sites, random_atom_types, find_atom_site_with_minimum_distance
from htsohm.pair_distance import min_pair_distance

def convert_positions_to_offset(x0, x1):
    offset = x1 - x0
    if offset > 0.5:
        return offset - 1.0
    elif offset < -0.5:
        return offset + 1.0
    else:
        return offset

def mutate_pos_to_new_pos_w_pbc(x0, x1, mutation_strength):
    """
    if given a position x0 and a second (usually random) position x1, then moves towards the nearest
    x1 (either within the current image, or across a PBC) a distance equal to the mutation_strength.
    """
    return (x0 + (convert_positions_to_offset(x0, x1)) * mutation_strength) % 1.

def move_sites(sitesl, uc_a, ms, distance, num_trials=100):
    for site in sitesl:
        other_sites = set(sitesl) - {site}
        good_pos = find_good_move_position(site, other_sites, uc_a, ms, distance, num_trials)
        if good_pos is not None:
            site.xyz = good_pos
        else:
            slog("Failed to move site")

def find_good_move_position(site, other_sites, uc_a, mutation_strength, distance, num_trials=100):
    """ attempts {num_trials} trials of finding a new move position for one site in a collection of
    sites, such that the minimum distance between points is maintained at > {distance}
    """
    for _ in range(num_trials):
        trial_pos = [mutate_pos_to_new_pos_w_pbc(r, random(), mutation_strength) for r in site.xyz]
        if min_pair_distance([s.xyz for s in other_sites] + [trial_pos]) > distance/uc_a:
            return trial_pos
    return None


def rebalance_charges(charges, charge_limits, exclude_indices=[]):
    """rebalances the list charges in-place so that the total charge equals 0.0 and the
    individual charges are withih the charge_limits.
    """

    delta_q = sum(charges)
    # reassign delta_q to other atoms
    other_charges = set(range(len(charges))) - set(exclude_indices)
    while not math.isclose(delta_q, 0.0, abs_tol=1e-12):
        i = choice(list(other_charges))
        i_new_q = min(max(charges[i] - delta_q, charge_limits[0]), charge_limits[1])
        delta_q += i_new_q - charges[i]
        charges[i] = i_new_q
        other_charges = other_charges - {i}

    return charges

def mutate_charge(charges, charge_index, mutation_strength, charge_limits):
    if len(charges) == 1:
        return [0.0]

    new_charges = charges.copy()
    new_q = perturb_unweighted(charges[charge_index], mutation_strength, charge_limits)
    new_charges[charge_index] = new_q
    return rebalance_charges(new_charges, charge_limits, exclude_indices=[charge_index])

def mutate_charges(charges, mutation_strength, charge_limits):
    new_charges = charges
    for i in range(len(new_charges)):
        new_charges = mutate_charge(new_charges, i, mutation_strength, charge_limits)

    return new_charges




def net_charge(atom_sites):
    return sum([e.q for e in atom_sites])

def perturb_unweighted(curr_val, mutation_strength, var_limits):
    # the / 2 is to make sure that if the max_change is the entire variable range (i.e. 100%
    # mutation strength), then if we are in the center of the range, this equivalent to generating
    # a new random number in the range.
    max_change = mutation_strength * (var_limits[1] - var_limits[0])
    new_val = curr_val + uniform(-max_change, max_change) / 2
    return min(max(new_val, var_limits[0]), var_limits[1])

def print_parent_child_diff(parent, child):
    slog("PARENT ID :\t{}".format(parent.id))
    slog("CHILD ID  :\t{}".format(child.id))
    slog("lattice constants: (%.2f, %.2f, %.2f) => (%.2f, %.2f, %.2f)" % (parent.a, parent.b, parent.c, child.a, child.b, child.c))
    slog("number of atoms: %.2f => %.2f" % (len(parent.atom_sites), len(child.atom_sites)))
    parent_ats = ", ".join(["(%.1f, %.1f)" % (ats.epsilon, ats.sigma) for ats in parent.atom_types])
    child_ats = ", ".join(["(%.1f, %.1f)" % (ats.epsilon, ats.sigma) for ats in child.atom_types])
    slog("atom types: %s => %s" % (parent_ats, child_ats))

def mutate_material(parent, config):
    child = parent.clone()
    child.parent = parent
    child.parent_id = parent.id

    perturb = set(config["perturb"])
    if config["perturb_type"] == "random":
        child.perturbation = choice(perturb)
        perturb = {child.perturbation}
    else:
        child.perturbation = "all"

    slog("Parent id: %d" % (child.parent_id))
    slog("Perturbing: %s [%s]" % (child.perturbation, perturb))
    ms = config["mutation_strength"]

    if config["number_of_atom_types"] > len(child.atom_types):
        num_atom_types_to_add = config["number_of_atom_types"] - len(child.atom_types)
        slog("Adding %d random atom types so we have number defined in the config" % num_atom_types_to_add)
        child.atom_types += random_atom_types(num_atom_types_to_add, config)

    if perturb & {"num_atoms"} and random() < ms:
        if random() < 0.5: # remove an atom
            if len(child.atom_sites) > config['num_atoms_limits'][0]:
                site_to_remove = choice(child.atom_sites)
                slog("Removing atom site: ", site_to_remove)
                removed_site = child.atom_sites.remove(site_to_remove)
                child.allq = rebalance_charges(child.allq, config['charge_limits'])
        else: # add an atom
            if len(child.atom_sites) < config['num_atoms_limits'][1]:
                slog("Adding atom site...")
                atom_position = find_atom_site_with_minimum_distance([s.xyz for s in child.atom_sites], config['minimum_site_distance'], child.a)
                if atom_position:
                    # note, the xyz coordinates from random_atom_sites will be replaced by the point found above
                    new_site = random_atom_sites(1, child.atom_types, q=[uniform(*config['charge_limits'])])[0]
                    new_site.xyz = atom_position
                    child.atom_sites.append(new_site)
                    child.allq = rebalance_charges(child.allq, config['charge_limits'], exclude_indices=[len(child.atom_sites) - 1])
                else:
                    slog("Failed to add a new atom.")

    if perturb & {"atom_type_assignments"}:
        for i, atom in enumerate(child.atom_sites):
            if random() < ms**2:
                new_atom_type = choice(child.atom_types)
                slog("Reassigning atom type for site %d from %d to %d" %
                    (i, child.atom_types.index(atom.atom_types), child.atom_types.index(new_atom_type)))
                atom.atom_types = new_atom_type

    if perturb & {"atom_types"}:
        sigl = config["sigma_limits"]
        epsl = config["epsilon_limits"]
        for at in child.atom_types:
            at.sigma = perturb_unweighted(at.sigma, ms, sigl)
            at.epsilon = perturb_unweighted(at.epsilon, ms, epsl)

    if perturb & {"lattice"}:
        ll = config["lattice_constant_limits"]
        trial_a = perturb_unweighted(child.a, ms, ll)
        child.a = max(trial_a, child.min_unit_cell_a(config['minimum_site_distance']))
        if config["lattice_cubic"]:
            child.b = child.a
            child.c = child.a
        else:
            child.b = perturb_unweighted(child.b, ms, ll)
            child.c = perturb_unweighted(child.c, ms, ll)

    if perturb & {"atom_sites"}:
        move_sites(child.atom_sites, child.a, ms, config['minimum_site_distance'])

    if perturb & {"charges"}:
        child.allq = mutate_charges(child.allq, ms, config['charge_limits'])

    # possibility that the material is unchanged

    print_parent_child_diff(parent, child)
    return child
