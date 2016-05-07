import os
import numpy as np

from htsohm.runDB_declarative import Base, RunData, session
from htsohm import simulate as sim

def count_bin(run_id, ml_bin, sa_bin, vf_bin):
    
    
    c0 = session.query(RunData).filter(RunData.run_id == run_id,              # Not dummy-tested
                                        RunData.methane_loading_bin == ml_bin,
                                        RunData.surface_area_bin == sa_bin,
                                        RunData.void_fraction_bin == vf_bin,
                                        RunData.dummy_test_result == None).count()
    c1 = session.query(RunData).filter(RunData.run_id == run_id,              # Passed dummy-test
                                        RunData.methane_loading_bin == ml_bin,
                                        RunData.surface_area_bin == sa_bin,
                                        RunData.void_fraction_bin == vf_bin,
                                        RunData.dummy_test_result == 'y').count()
#    c2 = session.query(RunData).filter(RunData.Run == run_id,              # Parent failed dummy-test, child not tested
#                                 RunData.Bin_ML == ml_bin,
#                                 RunData.Bin_SA == sa_bin,
#                                 RunData.Bin_VF == vf_bin,
#                                 RunData.D_pass == 'm').count()


    bin_count = c0 + c1 #+ c2

    return bin_count


def check_number_of_bins(run_id):

    wd = os.environ['HTSOHM_DIR']
    with open( wd + '/' + run_id + '.txt' ) as origin:
        for line in origin:
            if "Number of bins:" in line:
                bins = int( line.split()[3] )

    return bins


def count_all(run_id):

    bins = check_number_of_bins(run_id)
    all_counts = np.zeros([bins, bins, bins])

    for i in range(bins):
        for j in range(bins):
            for k in range(bins):
                
                b_count = count_bin(run_id, i, j, k)
                all_counts[i,j,k] = b_count

    return all_counts


def select_parents(run_id, children_per_generation, generation):

    bins = check_number_of_bins(run_id)
    counts = count_all(run_id)
    weights = np.zeros([bins, bins, bins])

    for i in range(bins):
        for j in range(bins):
            for k in range(bins):
                if counts[i,j,k] != 0.:
                    weights[i,j,k] = counts.sum() / counts[i,j,k]
    weights = weights / weights.sum()

    w_list = []
    id_list = []
    for i in range(bins):
        for j in range(bins):
            w_list = np.concatenate( [w_list, weights[i,j,:]] )
            for k in range(bins):
                bin_ids = []
                res = session.query(RunData).filter(RunData.run_id == run_id,
                                              RunData.methane_loading_bin == i,
                                              RunData.surface_area_bin == j,
                                              RunData.void_fraction_bin == k,
                                              RunData.dummy_test_result == None
                                              ).all()
                for item in res:
                    bin_ids.append(item.material_id)
                res = session.query(RunData).filter(RunData.run_id == run_id,
                                              RunData.methane_loading_bin == i,
                                              RunData.surface_area_bin == j,
                                              RunData.void_fraction_bin == k,
                                              RunData.dummy_test_result == 'y'
                                              ).all()
                for item in res:
                    bin_ids.append(item.id)
#                res = session.query(RunData).filter(RunData.Run == run_id,
#                                              RunData.Bin_ML == i,
#                                              RunData.Bin_SA == j,
#                                              RunData.Bin_VF == k,
#                                              RunData.D_pass == 'm').all()
#                for item in res:
#                    bin_ids.append(item.Mat)
                id_list = id_list + [bin_ids]

    first = generation * children_per_generation
    last = (generation + 1) * children_per_generation
    new_material_ids = np.arange(first, last)                   # IDs for next generation of materials
    new_material_primary_keys = []
    for i in new_material_ids:
        res = session.query(RunData).filter(RunData.run_id == run_id,
                                            RunData.material_id == i)
        for item in res:
            new_material_primary_keys.append(item.id)

    next_materials_list = []
    for i in new_material_primary_keys:

        parent_bin = np.random.choice(id_list, p=w_list)
        parent_id = np.random.choice(parent_bin)           # Select parent for new material
        next_material = [ i, parent_id ]
        next_materials_list.append( next_material )

    return next_materials_list


def add_parent_ids(run_id, next_materials_list):

    for i in next_materials_list:
        row = session.query(RunData).get(i[0])
        row.parent_id = i[1]
        session.commit()
