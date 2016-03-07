import os
import numpy as np

def find_missing(mat_dir):

    skip = len(mat_dir) + 1
    mats = [x[0][skip:] for x in os.walk(mat_dir)][1:]

    ml = np.genfromtxt(mat_dir + '/ch4_abs_cc_cc.txt', usecols=0, dtype=str)
    sa = np.genfromtxt(mat_dir + '/SAdata_m2_cc.txt', usecols=0, dtype=str)
    vf = np.genfromtxt(mat_dir + '/HVdata2col.txt', usecols=0, dtype=str)

    missing_ml = []
    missing_sa = []
    missing_vf = []

    for i in mats:
        if i not in ml:
            missing_ml.append(i)
        if i not in sa:
            missing_sa.append(i)
        if i not in vf:
            missing_vf.append(i)
		
        if missing_ml != []:
            print( 'Missing ' + len(missing_ml) + ' methane loading calculations.' )
    
        if missing_sa != []:
            print( 'Missing ' + len(missing_sa) + ' surface area calculations.' )
    
        if missing_vf != []:
            print( 'Missing ' + len(missing_vf) + ' void fraction calculations.' )
    
        if missing_ml == [] and missing_sa == [] and missing_vf == []:
            print('No missing data points!')

        return missing_ml, missing_sa, missing_vf
