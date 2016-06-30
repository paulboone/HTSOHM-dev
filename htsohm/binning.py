# related third party imports
import numpy as np
from sqlalchemy import func

# local application/library specific imports
from htsohm.runDB_declarative import Base, Material, session

def select_parents(run_id, children_per_generation, generation):
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
            func.count(Material.id), Material.methane_loading_bin, Material.surface_area_bin,
            Material.void_fraction_bin
        ) \
        .filter(Material.run_id == run_id, Material.dummy_test_result != 'fail') \
        .group_by(
            Material.methane_loading_bin, Material.surface_area_bin, Material.void_fraction_bin
        ).all()[1:]
    bins = [{"ML" : i[1], "SA" : i[2], "VF" : i[3]} for i in bins_and_counts]
    total = sum([i[0] for i in bins_and_counts])
    # ...then assigned a weight.
    weights = [i[0] / float(total) for i in bins_and_counts]

    ############################################################################
    # A parent-material is selected for each material in the next generation.
    next_generation = session \
        .query(Material) \
        .filter(Material.run_id == run_id, Material.generation == generation).all()

    for child in next_generation:
        # First, the bin is selected...
        parent_bin = np.random.choice(bins, p=weights)
        parent_query = session \
            .query(Material.id) \
            .filter(
                Material.run_id == run_id,
                Material.methane_loading_bin == parent_bin['ML'],
                Material.surface_area_bin == parent_bin['SA'],
                Material.void_fraction_bin == parent_bin['VF'],
                Material.dummy_test_result != 'fail'
            ).all()
        potential_parents = [i[0] for i in parent_query]
        # ...then a parent is select from the materials in that bin.
        parent_id = np.random.choice(potential_parents)
        child.parent_id = str(parent_id)

    return next_generation
