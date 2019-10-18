from random import choice, random, uniform, sample

from htsohm.generator.random import random_atom_sites, random_number_density

def random_position(x0, x1, mutation_strength):
    # get minimum distance between two points of (1) within the box, and (2) across the box boundary
    dx = min(abs(x0 - x1), 1 - abs(x0 - x1))
    if x0 > x1 and (x0 - x1) > 0.5: # then dx will be through boundary, so move right
        x2 = (x0 + mutation_strength * dx) % 1.
    if x0 >= x1 and (x0 - x1) < 0.5: # then dx will be in box, so move left
        x2 = x0 - mutation_strength * dx
    if x0 < x1 and (x1 - x0) > 0.5: # then dx will be through boundary, so move left
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
    print("PARENT UUID :\t{}".format(parent.uuid))
    print("CHILD UUID  :\t{}".format(child.uuid))
    print("lattice constants: (%.2f, %.2f, %.2f) => (%.2f, %.2f, %.2f)" % (ps.a, ps.b, ps.c, cs.a, cs.b, cs.c))
    print("number_density: %.2e => %.2e" % (parent.number_density, child.number_density))
    print("number of atoms: %.2f => %.2f" % (len(parent.structure.atom_sites),
                                             len(child.structure.atom_sites)))
    parent_ljs = ", ".join(["(%.1f, %.1f)" % (ljs.epsilon, ljs.sigma) for ljs in ps.lennard_jones])
    child_ljs = ", ".join(["(%.1f, %.1f)" % (ljs.epsilon, ljs.sigma) for ljs in cs.lennard_jones])
    print("lennard jones: %s => %s" % (parent_ljs, child_ljs))
    # print("FRAMEWORK NET CHARGE :\t{}".format(sum([e.q for e in cs.atom_sites])))

def mutate_material(parent, config):
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
    ms = config["mutation_strength"]

    if perturb & {"atom_types"}:
        sigl = config["sigma_limits"]
        epsl = config["epsilon_limits"]
        for at in cs.lennard_jones:
            at.sigma = perturb_unweighted(at.sigma, ms * (sigl[1] - sigl[0]), sigl)
            at.epsilon = perturb_unweighted(at.epsilon, ms * (epsl[1] - epsl[0]), epsl)

    if perturb & {"lattice"}:
        ll = config["lattice_constant_limits"]
        cs.a = perturb_unweighted(cs.a, ms * (ll[1] - ll[0]), ll)
        cs.b = perturb_unweighted(cs.b, ms * (ll[1] - ll[0]), ll)
        cs.c = perturb_unweighted(cs.c, ms * (ll[1] - ll[0]), ll)
        child.number_density = len(cs.atom_sites) / cs.volume

    if perturb & {"density"}:
        ndl = config["number_density_limits"]
        child.number_density = perturb_unweighted(child.number_density, (ndl[1] - ndl[0])*ms, ndl)

    if perturb & {"atom_sites"}:
        for a in cs.atom_sites:
            a.x = random_position(a.x, random(), ms)
            a.y = random_position(a.y, random(), ms)
            a.z = random_position(a.z, random(), ms)

    # add / remove atoms if density has changed
    if "fix_atoms" in config:
        number_of_atoms = config['fix_atoms']
    else:
        number_of_atoms = max(1, round(child.number_density * child.structure.volume))

    if number_of_atoms < len(cs.atom_sites):
        cs.atom_sites = sample(cs.atom_sites, number_of_atoms)
    elif number_of_atoms > len(cs.atom_sites):
        cs.atom_sites += random_atom_sites(number_of_atoms - len(cs.atom_sites), cs.lennard_jones)

    print_parent_child_diff(parent, child)
    return child
