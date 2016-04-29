# /usr/bin/env python
import sys
# import os
# parentdir = os.path.dirname(__file__)
# print(__file__)
# print(parentdir)
# exit()
#
# sys.path.insert(0,os.path.dirname(__file__))
#
# sys.path.insert(os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir))
# sys.path.insert(0,'')

from htsohm.HTSOHM import HTSOHM

if __name__ == "__main__":

    HTSOHM(int(sys.argv[1]),
           int(sys.argv[2]),
           float(sys.argv[3]),
           int(sys.argv[4]))
