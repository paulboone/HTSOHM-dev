from math import sqrt

from sqlalchemy.sql import func, or_
from sqlalchemy.orm.exc import FlushError

from htsohm import config
from htsohm.db import session, Material, MutationStrength
from htsohm.simulate import run_all_simulations
from htsohm.material_files import write_seed_definition_files, write_child_definition_files

def materials_in_generation(run_id, generation):
    return session.query(Material).filter(
        Material.run_id == run_id,
        Material.generation == generation
    ).count()

def last_generation(run_id):
    return session.query(func.max(Material.generation)).filter(
        Material.run_id == run_id,
    )[0][0]

def evaluate_convergence(run_id):
    '''Counts number of materials in each bin and returns variance of these counts.'''
    bin_counts = session \
        .query(func.count(Material.id)) \
        .filter(Material.run_id == run_id) \
        .group_by(
            Material.methane_loading_bin, Material.surface_area_bin, Material.void_fraction_bin
        ).all()
    bin_counts = [i[0] for i in bin_counts]    # convert SQLAlchemy result to list
    variance = sqrt( sum([(i - (sum(bin_counts) / len(bin_counts)))**2 for i in bin_counts]) / len(bin_counts))
    return variance


def select_parent(run_id, max_generation, generation_limit):
    """Use bin-counts to preferentially select a list of rare parents.

    Each bin contains some number of materials, and those bins with the fewers materials represent
    the most rare structure-property combinations. These rare materials are preferred as parents
    for new materials, because their children are most likely to display unique properties. This
    function first calculates a `weight` for each bin, based on the number of constituent
    materials. These weights affect the probability of selecting a parent from that bin. Once a bin
    is selected, a parent is randomly-selected from those materials within that bin.
    """

    # Each bin is counted...
    bins_and_counts = session \
        .query(
            func.count(Material.id),
            Material.methane_loading_bin,
            Material.surface_area_bin,
            Material.void_fraction_bin
        ) \
        .filter(
            Material.run_id == run_id,
            or_(Material.retest_passed == True, Material.retest_passed == None),
            Material.generation <= max_generation,
            Material.generation_index < generation_limit,
        ) \
        .group_by(
            Material.methane_loading_bin, Material.surface_area_bin, Material.void_fraction_bin
        ).all()[1:]
    bins = [{"ML" : i[1], "SA" : i[2], "VF" : i[3]} for i in bins_and_counts]
    total = sum([i[0] for i in bins_and_counts])
    # ...then assigned a weight.
    weights = [i[0] / float(total) for i in bins_and_counts]

    parent_bin = np.random.choice(bins, p=weights)
    parent_query = session \
        .query(Material.id) \
        .filter(
            Material.run_id == run_id,
            or_(Material.retest_passed == True, Material.retest_passed == None),
            Material.methane_loading_bin == parent_bin["ML"],
            Material.surface_area_bin == parent_bin["SA"],
            Material.void_fraction_bin == parent_bin["VF"],
            Material.generation <= max_generation,
            Material.generation_index < generation_limit,
        ).all()
    potential_parents = [i[0] for i in parent_query]

    return int(np.random.choice(potential_parents))


def mutate(run_id, generation, parent):
    """Retrieve the latest mutation_strength for the parent, or calculate it if missing.

    In the event that a particular bin contains parents whose children exhibit radically
    divergent properties, the strength parameter for the bin is modified. In order to determine
    which bins to adjust, the script refers to the distribution of children in the previous
    generation which share a common parent. The criteria follows:
     ________________________________________________________________
     - if none of the children share  |  halve strength parameter
       the parent's bin               |
     - if the fraction of children in |
       the parent bin is < 10%        |
     _________________________________|_____________________________
     - if the fraction of children in |  double strength parameter
       the parent bin is > 50%        |
     _________________________________|_____________________________
    """

    mutation_strength_key = [run_id, generation] + parent.bin
    mutation_strength = session.query(MutationStrength).get(mutation_strength_key)

    if mutation_strength:
        print("Mutation strength already calculated for this bin and generation.")
    else:
        print("Calculating mutation strength...")
        mutation_strength = MutationStrength.get_prior(*mutation_strength_key).clone()
        mutation_strength.generation = generation

        try:
            fraction_in_parent_bin = parent.calculate_percent_children_in_bin()
            if fraction_in_parent_bin < 0.1:
                mutation_strength.strength *= 0.5
            elif fraction_in_parent_bin > 0.5 and mutation_strength.strength <= 0.5:
                mutation_strength.strength *= 2
        except ZeroDivisionError:
            print("No prior generation materials in this bin with children.")

        try:
            session.add(mutation_strength)
            session.commit()
        except FlushError as e:
            print("Somebody beat us to saving a row with this generation. That's ok!")
            # it's ok b/c this calculation should always yield the exact same result!

    return mutation_strength.strength

def retest(m_orig, retests, tolerance):
    """Recalculate material structure-properties to prevent statistical errors.

    Because methane loading, surface area, and helium void fractions are calculated using
    statistical methods (namely grand canonic Monte Carlo simulations) they are susceptible
    to statistical errors. To mitigate this, after a material has been selected as a potential
    parent, it's combination of structure-properties is resimulated some number of times and
    compared to the initally-calculated material. If the resimulated values differ from the
    initially-calculated value beyond an accpetable tolerance, the material fails the `dummy-test`
    and is flagged, preventing it from being used to generate new materials in the future.
    """

    m = m_orig.clone()
    run_all_simulations(m)

    # requery row from database, in case someone else has changed it, and lock it
    # if the row is presently locked, this method blocks until the row lock is released
    session.refresh(m_orig, lockmode='update')
    if m_orig.retest_num < retests:
        m_orig.retest_methane_loading_sum += m.ml_absolute_volumetric_loading
        m_orig.retest_surface_area_sum += m.sa_volumetric_surface_area
        m_orig.retest_void_fraction_sum += m.vf_helium_void_fraction
        m_orig.retest_num += 1

        if m_orig.retest_num == retests:
            m_orig.retest_passed = m.calculate_retest_result(tolerance)
    else:
        pass
        # otherwise our test is extra / redundant and we don't save it

    session.commit()


def worker_run_loop(run_id):
    gen = last_generation(run_id) or 0

    converged = False
    while not converged:
        size_of_generation = config['children-per-generation']

        while materials_in_generation(run_id, gen) < size_of_generation:
            if gen == 0:
                print("writing new seed...")
                material = write_seed_definition_files(run_id, config['number-of-atom-types'])
            else:
                print("selecting a parent / running retests on parent / mutating / simulating")
                parent_id = select_parent(run_id, max_generation=(gen - 1),
                                                  generation_limit=config['children-per-generation'])

                parent = session.query(Material).get(parent_id)

                # run retests until we've run enough
                while parent.retest_passed is None:
                    print("running retest...")
                    retest(parent, config['retests']['number'], config['retests']['tolerance'])
                    session.refresh(parent)

                if not parent.retest_passed:
                    print("parent failed retest. restarting with parent selection.")
                    continue

                mutation_strength = mutate(run_id, gen, parent)
                material = write_child_definition_files(run_id, parent_id, gen, mutation_strength)

            run_all_simulations(material)
            session.add(material)
            session.commit()

            material.generation_index = material.calculate_generation_index()
            if material.generation_index < config['children-per-generation']:
                session.add(material)
            else:
                # delete excess rows
                session.delete(material)
            session.commit()
        gen += 1

        # no convergance test at present!
