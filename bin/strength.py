import numpy as np
import os
import sys
#from binning import *

#HTSOHM_dir
#p_dir
#bins

def calcS(HTSOHM_dir, p_dir, bins):

    sys.path.insert(0, HTSOHM_dir)

    from binning import *

    gp_dir = 'gen' + str(int(p_dir[4:]) - 1)

    

