# standard library imports
import os

# related third party imports
import numpy as np
import yaml

# local application/library specific imports
from htsohm.runDB_declarative import Base, Material, session
from htsohm.utilities import read_config_file

def count_bin(run_id, ml_bin, sa_bin, vf_bin):
    """Returns number of materials in a particular bin."""
    bin_count = session.query(Material).filter(Material.run_id == run_id, Material.methane_loading_bin == ml_bin,
        Material.surface_area_bin == sa_bin, Material.void_fraction_bin == vf_bin,
        Material.dummy_test_result.in_(('none', 'pass'))).count()
    return bin_count

def count_all(run_id):
    """Returns an array containing the number of materials in every bin."""
    config = read_config_file(run_id)
    bins = config["number-of-bins"]
    all_counts = np.zeros([bins, bins, bins])
    for i in range(bins):
        for j in range(bins):
            for k in range(bins):
                bin_count = count_bin(run_id, i, j, k)
                all_counts[i,j,k] = bin_count
    return all_counts

def select_parents(run_id, children_per_generation, generation):
    """Use bin-counts to preferentially select a list of rare parents.

    Each bin contains some number of materials, and those bins with the fewers materials represent
    the most rare structure-property combinations. These rare materials are preferred as parents
    for new materials, because their children are most likely to display unique properties. This
    function first calculates a `weight` for each bin, based on the number of constituent
    materials. These weights affect the probability of selecting a parent from that bin. Once a bin
    is selected, a parent is randomly-selected from those materials within that bin."""

    # Each bin is counted, then assigned a weight
    config = read_config_file(run_id)
    bins = config["number-of-bins"]
    counts = count_all(run_id)
    weights = np.zeros([bins, bins, bins])
    for i in range(bins):
        for j in range(bins):
            for k in range(bins):
                if counts[i,j,k] != 0.:
                    weights[i,j,k] = counts.sum() / counts[i,j,k]
    weights = weights / weights.sum()

    ############################################################################
    # The 3-dimensional list of weights is converted into a 1-dimensional list,
    # and another 1-dimensional list of material-IDs is also generated.
    weight_list = []
    id_list = []
    for i in range(bins):
        for j in range(bins):
            weight_list = np.concatenate([weight_list, weights[i,j,:]])
            for k in range(bins):
                bin_ids = []
                materials = session.query(Material).filter(Material.run_id == run_id, Material.methane_loading_bin == i,
                    Material.surface_area_bin == j, Material.void_fraction_bin == k,
                    Material.dummy_test_result.in_(['none', 'pass'])).all()
                for material in materials:
                    bin_ids.append(material.id)
                id_list = id_list + [bin_ids]

    ############################################################################
    # A parent-material is selected for each material in the next generation.
    next_generation = session.query(Material).filter(Material.run_id == run_id, Material.generation == generation).all()

    for material in next_generation:
        parents_list = np.random.choice(id_list, p=weight_list)  # weighted sampling function
        parent_id = np.random.choice(parents_list)
        material.parent_id = str(parent_id)

    return next_generation
