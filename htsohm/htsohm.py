import os
import sys
from math import sqrt
from datetime import datetime

import numpy as np
from sqlalchemy.sql import func, or_
from sqlalchemy.orm.exc import FlushError
import sqlalchemy.exc
import yaml

import htsohm
from htsohm import config
from htsohm.db import session, Material, MutationStrength
from htsohm.material_files import generate_pseudo_material, mutate_pseudo_material
from htsohm import simulation

def materials_in_generation(run_id, generation):
    """Count number of materials in a generation.

    Args:
        run_id (str): identification string for run.
        generation (int): iteration in overall bin-mutate-simulate rountine.

    Returns:
        Number(int) of materials in a particular generation that are present in
        the database (the final step in bin-mutate-simulate routine).

    """
    return session.query(Material).filter(
        Material.run_id == run_id,
        Material.generation == generation
    ).count()

def last_generation(run_id):
    """Finds latest generation present in database.

    Args:
        run_id (str): identification string for run.

    Returns:
        Last generation(int) to be included in database.

    """
    return session.query(func.max(Material.generation)).filter(
        Material.run_id == run_id,
    )[0][0]

def calc_bin(value, bound_min, bound_max, bins):
    """Find bin in parameter range.

    Args:
        value (float): some value, the result of a simulation.
        bound_min (float): lower limit, defining the parameter-space.
        bound_max (float): upper limit, defining the parameter-space.
        bins (int): number of bins used to subdivide parameter-space.

    Returns:
        Bin(int) corresponding to the input-value.

    """
    step = (bound_max - bound_min) / bins
    assigned_bin = (value - bound_min) // step
    assigned_bin = min(assigned_bin, bins-1)
    assigned_bin = max(assigned_bin, 0)
    return int(assigned_bin)

def select_parent(run_id, max_generation, generation_limit):
    """Use bin-counts to preferentially select a list of 'rare' parents.

    Args:
        run_id (str): identification string for run.
        max_generation (int): latest generation to include when counting number
            of materials ub each bin.
        generation_limit (int): number of materials to query in each generation
            (as materials are added to database they are assigned an index
            within the generation to bound the number of materials in each
            generation).

    Returns:
        The material id(int) corresponding to some parent-material selected
        from database with a bias favoring materials in bins with the lowest
        counts.

    """
    simulations = config['material_properties']
    queries = []
    if 'gas_adsorption' in simulations:
        queries.append( getattr(Material, 'gas_adsorption_bin') )
    if 'surface_area' in simulations:
        queries.append( getattr(Material, 'surface_area_bin') )
    if 'helium_void_fraction' in simulations:
        queries.append( getattr(Material, 'void_fraction_bin') )

    # Each bin is counted...
    bins_and_counts = session \
        .query(
            func.count(Material.id),
            *queries
        ) \
        .filter(
            Material.run_id == run_id,
            or_(Material.retest_passed == True, Material.retest_passed == None),
            Material.generation <= max_generation,
            Material.generation_index < generation_limit,
        ) \
        .group_by(*queries).all()[1:]

    bins = []
    for i in bins_and_counts:
        bin = {}
        for j in range(len(queries)):
            bin[queries[j]] = i[j + 1]
        bins.append(bin)
    total = sum([i[0] for i in bins_and_counts])
    # ...then assigned a weight.
    weights = [ total / float(i[0]) for i in bins_and_counts ]
    normalized_weights = [ weight / sum(weights) for weight in weights ]
    parent_bin = np.random.choice(bins, p = normalized_weights)
    
    parent_queries = [i == parent_bin[i] for i in queries]
    parent_query = session \
        .query(Material.id) \
        .filter(
            Material.run_id == run_id,
            or_(Material.retest_passed == True, Material.retest_passed == None),
            *parent_queries,
            Material.generation <= max_generation,
            Material.generation_index < generation_limit,
        ).all()
    potential_parents = [i[0] for i in parent_query]
    return int(np.random.choice(potential_parents))

def run_all_simulations(material, pseudo_material):
    """Simulate helium void fraction, gas loading, and surface area.

    Args:
        material (sqlalchemy.orm.query.Query): material to be analyzed.

    Depending on properties specified in config, adds simulated data for helium
    void fraction, gas loading, heat of adsorption, surface area, and
    corresponding bins to row in database corresponding to the input-material.
        
    """
    simulations = config['material_properties']

    ############################################################################
    # run helium void fraction simulation
    if 'helium_void_fraction' in simulations:
        results = simulation.helium_void_fraction.run(
            material.run_id, pseudo_material)
        material.update_from_dict(results)
        material.void_fraction_bin = calc_bin(
            material.vf_helium_void_fraction,
            *config['helium_void_fraction']['limits'],
            config['number_of_convergence_bins']
        )
    else:
        material.void_fraction_bin = 0
    ############################################################################
    # run gas loading simulation
    if 'gas_adsorption_0' in simulations and 'gas_adsorption_1' not in simulations:
        arguments = [material.run_id, pseudo_material]
        if 'helium_void_fraction' in simulations:
            arguments.append(material.vf_helium_void_fraction)
        results = simulation.gas_adsorption_0.run(*arguments)
        material.update_from_dict(results)
        material.gas_adsorption_bin = calc_bin(
            material.ga0_absolute_volumetric_loading,
            *config['gas_adsorption_0']['limits'],
            config['number_of_convergence_bins']
        )
    elif 'gas_adsorption_0' in simulations and 'gas_adsorption_1' in simulations:
        arguments = [material.run_id, pseudo_material]
        if 'helium_void_fraction' in simulations:
            arguments.append(material.vf_helium_void_fraction)
        results = simulation.gas_adsorption_0.run(*arguments)
        material.update_from_dict(results)
        results = simulation.gas_adsorption_1.run(*arguments)
        material.update_from_dict(results)
        material.gas_adsorption_bin = calc_bin(
            abs(material.ga0_absolute_volumetric_loading - material.ga1_absolute_volumetric_loading),
            *config['gas_adsorption_0']['limits'],
            config['number_of_convergence_bins']
        )
    elif 'gas_adsorption_0' not in simulations:
        material.gas_adsorption_bin = 0
    ############################################################################
    # run surface area simulation
    if 'surface_area' in simulations:
        results = simulation.surface_area.run(
                material.run_id, pseudo_material)
        material.update_from_dict(results)
        material.surface_area_bin = calc_bin(
            material.sa_volumetric_surface_area,
            *config['surface_area']['limits'],
            config['number_of_convergence_bins']
        )
    else:
        material.surface_area_bin = 0

def retest(m_orig, retests, tolerance, pseudo_material):
    """Reproduce simulations  to prevent statistical errors.

    Args:
        m_orig (sqlalchemy.orm.query.Query): material to retest.
        retests (int): number of times to reproduce each simulation.
        tolerance (float): acceptance criteria as percent deviation from
            originally calculated value.

    Queries database to determine if there are remained retests to  be run.
    Updates row in database with total number of retests and results.

    """
    m = m_orig.clone()
    run_all_simulations(m, pseudo_material)
    print('\n\nRETEST_NUM :\t%s' % m_orig.retest_num)
    print('retests :\t%s' % retests)

    simulations = config['material_properties']

    # requery row from database, in case someone else has changed it, and lock it
    # if the row is presently locked, this method blocks until the row lock is released
    session.refresh(m_orig, lockmode='update')
    if m_orig.retest_num < retests:
        if 'gas_adsorption_0' in simulations:
            m_orig.retest_gas_adsorption_0_sum += m.ga0_absolute_volumetric_loading
        if 'gas_adsorption_1' in simulations:
            m_orig.retest_gas_adsorption_1_sum += m.ga1_absolute_volumetric_loading
        if 'surface_area' in simulations:
            m_orig.retest_surface_area_sum += m.sa_volumetric_surface_area
        if 'helium_void_fraction' in simulations:
            m_orig.retest_void_fraction_sum += m.vf_helium_void_fraction
        m_orig.retest_num += 1
        session.commit()        

        if m_orig.retest_num == retests:
            try:
                m_orig.retest_passed = m.calculate_retest_result(tolerance)
                print('\nRETEST_PASSED :\t%s' % m_orig.retest_passed)
            except ZeroDivisionError as e:
                print('WARNING: ZeroDivisionError - material.calculate_retest_result(tolerance)')

    else:
        pass
        # otherwise our test is extra / redundant and we don't save it

    session.commit()

def get_mutation_strength(run_id, generation, parent):
    """Query mutation_strength for bin and adjust as necessary.

    Args:
        run_id (str): identification string for run.
        generation (int): iteration in bin-mutate-simulate routine.
        parent (sqlalchemy.orm.query.Query): parent-material corresponding to
            the bin being queried.

    Returns:
        mutation_strength.strength (float): mutation strength to be used for
        parents in the bin being queried. If the fraction of children from
        previous generation which populate the SAME bin as their parent is LESS
        THAN 10% then the mutation strength is REDUCED BY HALF. If the fraction
        of these children populating the SAME bin as their parent is GREATER
        THAN 50% then the mutation strength is INCREASED BY 200%.

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
        except (FlushError, sqlalchemy.exc.IntegrityError) as e:
            print("Somebody beat us to saving a row with this generation. That's ok!")
            session.rollback()
            # it's ok b/c this calculation should always yield the exact same result!
    sys.stdout.flush()
    return mutation_strength.strength

def evaluate_convergence(run_id, generation):
    '''Determines convergence by calculating variance of bin-counts.
    
    Args:
        run_id (str): identification string for run.
        generation (int): iteration in bin-mutate-simulate routine.

    Returns:
        bool: True if variance is less than or equal to cutt-off criteria (so
            method will continue running).
    '''
    simulations = config['material_properties']
    query_group = []
    if 'gas_adsorption' in simulations:
        query_group.append( getattr(Material, 'gas_adsorption_bin') )
    if 'surface_area' in simulations:
        query_group.append( getattr(Material, 'surface_area_bin') )
    if 'helium_void_fraction' in simulations:
        query_group.append( getattr(Material, 'void_fraction_bin') )

    bin_counts = session \
        .query(func.count(Material.id)) \
        .filter(
            Material.run_id == run_id, Material.generation < generation,
            Material.generation_index < config['children_per_generation']
        ) \
        .group_by(*query_group).all()
    bin_counts = [i[0] for i in bin_counts]    # convert SQLAlchemy result to list
    variance = sqrt( sum([(i - (sum(bin_counts) / len(bin_counts)))**2 for i in bin_counts]) / len(bin_counts))
    print('\nCONVERGENCE:\t%s\n' % variance)
    sys.stdout.flush()
    return variance <= config['convergence_cutoff_criteria']

def worker_run_loop(run_id):
    """
    Args:
        run_id (str): identification string for run.

    Writes seed generation and simulates properties, then manages overall
    bin-mutate-simualte routine until convergence cutt-off or maximum
    number of generations is reached.

    """
    gen = last_generation(run_id) or 0

    converged = False
    while not converged:
        print(
                (
                    '=======================================================\n'
                    'GENERATION {0}\n'
                    '=======================================================\n'
                ).format(gen)
            )
        size_of_generation = config['children_per_generation']

        while materials_in_generation(run_id, gen) < size_of_generation:
            if gen == 0:
                print("writing new seed...")
                material, pseudo_material = generate_pseudo_material(
                        run_id, config['number_of_atom_types'])
                pseudo_material.dump()
            else:
                print("selecting a parent / running retests on parent / mutating / simulating")
                parent_id = select_parent(run_id, max_generation=(gen - 1),
                                                  generation_limit=config['children_per_generation'])

                parent_material = session.query(Material).get(parent_id)

                # run retests until we've run enough
                while parent_material.retest_passed is None:
                    print("running retest...")
                    print("Date :\t%s" % datetime.now().date().isoformat())
                    print("Time :\t%s" % datetime.now().time().isoformat())
                    parent_pseudo_material_path = os.path.join(
                            os.path.dirname(os.path.dirname(htsohm.__file__)),
                            run_id, 'pseudo_materials', '{0}.yaml'.format(
                                parent_material.uuid)
                        )
                    with open(parent_pseudo_material_path) as parent_file:
                        parent_pseudo_material = yaml.load(parent_file)
                    retest(
                            parent_material, config['retests']['number'],
                            config['retests']['tolerance'], parent_pseudo_material
                        )
                    session.refresh(parent_material)

                if not parent_material.retest_passed:
                    print("parent failed retest. restarting with parent selection.")
                    continue

                mutation_strength = get_mutation_strength(run_id, gen, parent_material)
                material, pseudo_material = mutate_pseudo_material(
                        parent_material, parent_pseudo_material, mutation_strength, gen)
                pseudo_material.dump()
            run_all_simulations(material, pseudo_material)
            session.add(material)
            session.commit()

            material.generation_index = material.calculate_generation_index()
            if material.generation_index < config['children_per_generation']:
                session.add(material)
            else:
                # delete excess rows
                # session.delete(material)
                pass
            session.commit()
            sys.stdout.flush()
        gen += 1
        converged = evaluate_convergence(run_id, gen)
