#!/usr/bin/env python3

import sys

from htsohm.htsohm_hpc import start_run

if __name__ == "__main__":
    start_run(int(sys.argv[1]),
              int(sys.argv[2]),
              float(sys.argv[3]),
              int(sys.argv[4]))
