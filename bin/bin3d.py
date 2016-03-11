import numpy as np


def bin3d(p_dir, bins):

    ml_file = p_dir + '/ch4_abs_cc_cc.txt'
    sa_file = p_dir + '/SAdata_m2_cc.txt'
    vf_file = p_dir + '/HVdata2col.txt'

    # bins data in one dimension (ml, sa, or vf)
    def bin_data(data_file, e_min, e_max, bins):

        name = np.genfromtxt(data_file, usecols=0, dtype=str)
        valu = np.genfromtxt(data_file, usecols=1, dtype=float)

        step = e_max / bins
        edge = np.arange(e_min, e_max + step, step)

        IDs = [[] for i in range(bins)]

        for i in range(bins):

            lower = edge[i]
            upper = edge[i + 1]

            for j in range(len(valu)):
                if valu[j] >= lower and valu[j] < upper:
                    IDs[i].append(name[j])

        return IDs

    ml_IDs = bin_data(ml_file, 0, 400., bins)
    sa_IDs = bin_data(sa_file, 0, 4500., bins)
    vf_IDs = bin_data(vf_file, 0, 1., bins)

    bin_IDs = [[[[] for i in range(bins)] for j in range(
        bins)] for k in range(bins)]  # make empty 3D list

    freq = np.zeros([bins, bins, bins])

    for i in range(bins):
        for j in range(bins):
            for k in range(bins):
                for l in ml_IDs[i]:
                    for m in sa_IDs[j]:
                        if l == m:
                            for n in vf_IDs[k]:
                                if m == n:
                                    bin_IDs[i][j][k].append(n)

    # 3D array containing number of materials per bin
    bin_counts = np.empty([bins, bins, bins])

    for i in range(bins):
        for j in range(bins):
            for k in range(bins):
                bin_counts[i, j, k] = len(bin_IDs[i][j][k])

    return bin_counts, bin_IDs
