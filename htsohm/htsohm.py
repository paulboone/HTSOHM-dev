
from sqlalchemy.sql import func
from sqlalchemy.orm.exc import FlushError

from htsohm import config
from htsohm.binning import select_parent
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
        m_orig.retest_methane_loading_sum += m.absolute_volumetric_loading
        m_orig.retest_surface_area_sum += m.volumetric_surface_area
        m_orig.retest_void_fraction_sum += m.helium_void_fraction
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
