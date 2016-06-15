#!/usr/bin/env python3

import sys

from htsohm.htsohm import htsohm

if __name__ == "__main__":
    htsohm(int(sys.argv[1]),
           int(sys.argv[2]),
           float(sys.argv[3]),
           int(sys.argv[4]))
