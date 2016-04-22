import numpy as np

from htsohm.runDB_declarative import Base, RunData
from htsohm.simulate import CreateSession, GetValue, AddRows, UpdateTable


def CountBin(run_ID, ML_bin, SA_bin, VF_bin):
    
    s = CreateSession()
    c0 = s.query(RunData).filter(RunData.run_id == run_ID,              # Not dummy-tested
                                 RunData.methane_loading_bin == ML_bin,
                                 RunData.surface_area_bin == SA_bin,
                                 RunData.void_fraction_bin == VF_bin,
                                 RunData.dummy_test_result == None).count()
    c1 = s.query(RunData).filter(RunData.run_id == run_ID,              # Passed dummy-test
                                 RunData.methane_loading_bin == ML_bin,
                                 RunData.surface_area_bin == SA_bin,
                                 RunData.void_fraction_bin == VF_bin,
                                 RunData.dummy_test_result == 'y').count()
#    c2 = s.query(RunData).filter(RunData.Run == run_ID,              # Parent failed dummy-test, child not tested
#                                 RunData.Bin_ML == ML_bin,
#                                 RunData.Bin_SA == SA_bin,
#                                 RunData.Bin_VF == VF_bin,
#                                 RunData.D_pass == 'm').count()


    BinCount = c0 + c1 #+ c2

    return BinCount


def CountAll(run_ID):

    bins = int(run_ID[-1])
    AllCounts = np.zeros([bins, bins, bins])

    for i in range(bins):
        for j in range(bins):
            for k in range(bins):
                
                b_count = CountBin(run_ID, i, j, k)
                AllCounts[i,j,k] = b_count

    return AllCounts


def SelectParents(run_ID, children_per_generation, generation):

    s = CreateSession()

    bins = int(run_ID[-1])
    counts = CountAll(run_ID)
    weights = np.zeros([bins, bins, bins])

    for i in range(bins):
        for j in range(bins):
            for k in range(bins):    
                if counts[i,j,k] != 0.:
                    weights[i,j,k] = counts.sum() / counts[i,j,k]
    weights = weights / weights.sum()

    w_list = []
    ID_list = []
    for i in range(bins):
        for j in range(bins):
            w_list = np.concatenate( [w_list, weights[i,j,:]] )
            for k in range(bins):
                bin_IDs = []
                res = s.query(RunData).filter(RunData.run_id == run_ID,
                                              RunData.methane_loading_bin == i,
                                              RunData.surface_area_bin == j,
                                              RunData.void_fraction_bin == k,
                                              RunData.dummy_test_result == None
                                              ).all()
                for item in res:
                    bin_IDs.append(item.material_id)
                res = s.query(RunData).filter(RunData.run_id == run_ID,
                                              RunData.methane_loading_bin == i,
                                              RunData.surface_area_bin == j,
                                              RunData.void_fraction_bin == k,
                                              RunData.dummy_test_result == 'y'
                                              ).all()
                for item in res:
                    bin_IDs.append(item.material_id)
#                res = s.query(RunData).filter(RunData.Run == run_ID,
#                                              RunData.Bin_ML == i,
#                                              RunData.Bin_SA == j,
#                                              RunData.Bin_VF == k,
#                                              RunData.D_pass == 'm').all()
#                for item in res:
#                    bin_IDs.append(item.Mat)
                ID_list = ID_list + [bin_IDs]

    first = generation * children_per_generation
    last = (generation + 1) * children_per_generation
    new_mat_IDs = np.arange(first, last)                   # IDs for next generation of materials

    for i in new_mat_IDs:

        p_bin = np.random.choice(ID_list, p=w_list)
        p_ID = np.random.choice(p_bin)                     # Select parent for new material

        _id =GetValue(run_ID, p_ID, "id")

        AddRows(run_ID, [i])
        data = {'parent_id': _id}
        UpdateTable(run_ID, i, data)
   

#def CheckConvergance(run_ID):
#
#    counts = CountAll(run_ID)
       
