# /usr/bin/env python

import sys
import os
import datetime 

import numpy as np

import generate as gen
import simulate as sim
import binning as bng
import mutate as mut


def HTSOHM(children_per_generation,    # number of materials per generation
           number_of_atomtypes,        # number of atom-types per material
           strength_0,                 # intial strength parameter
           number_of_bins,             # number of bins for analysis
           max_generations=20):        # maximum number of generations

    # Start run (see DD.MM.YYYY_HH.MM.SS_CpG.NoA_S0_NoB.txt for parameters)
    sys.path.insert(0, os.environ['SRC_DIR'])

    start = datetime.datetime.now()
    run_ID = ( "%s.%s.%s_%s.%s.%s_%s.%s_%s_%s" %
               (start.day, start.month, start.year,
                start.hour, start.minute, start.second,
                children_per_generation, number_of_atomtypes,
                strength_0,
                number_of_bins))

    wd = os.environ['HTSOHM_DIR']      # specify working directory          

    run_file = open( wd + '/' + run_ID + '.txt', "w")
    run_file.write( "Date:\t\t\t\t%s:%s:%s\n" % (start.day, start.month,
                                                 start.year) +
                    "Time:\t\t\t\t%s:%s:%s\n" % (start.hour, start.minute, 
                                                 start.second) +
                    "Children per generation:\t%s\n" % (
                                                 children_per_generation) +
                    "Number of atom-types:\t\t%s\n" % (number_of_atomtypes) +
                    "Initial mutation strength:\t%s\n" % (strength_0) +
                    "Number of bins:\t\t\t%s\n" % (number_of_bins))
    run_file.close()


    # SEED (GENERATION = 0)
    generation = 0
    gen.generate(children_per_generation, number_of_atomtypes, run_ID)
    sim.AddRows( run_ID, GenIDs(generation, children_per_generation) )
    sim.simulate( run_ID, GenIDs(generation, children_per_generation) )


    # FIRST GENERATION
    generation = 1                     # `Generation` counter
    # Select parents, add IDs to database...
    bng.SelectParents(run_ID, children_per_generation, generation)
    sim.DummyTest(run_ID, generation)
    
    mut.FirstS(run_ID, strength_0)     # Create strength-parameter array `run_ID`.npy
    mut.mutate(run_ID, generation)     # Create first generation of child-materials

    sim.AddRows( run_ID, GenIDs(generation, children_per_generation) )
    sim.simulate( run_ID, GenIDs(generation, children_per_generation) )


    # SECOND GENERATION, AND ON...
    NextGens = np.arange(2, max_generations)
    for i in NextGens:
        generation = i

        bng.SelectParents(run_ID, children_per_generation, generation)
        sim.DummyTest(run_ID, generation)
        mut.CalculateS(run_ID, generation)
        mut.mutate(run_ID, generation)
        sim.AddRows( run_ID, GenIDs(generation, children_per_generation) )
        sim.simulate( run_ID, GenIDs(generation, children_per_generation) )

def GenIDs(generation, children_per_generation):
    first = generation * children_per_generation
    last = (generation + 1) * children_per_generation
    GenIDs = np.arange(first, last)

    return GenIDs


if __name__ == "__main__":
    import sys
    HTSOHM(int(sys.argv[1]),
           int(sys.argv[2]),
           float(sys.argv[3]),
           int(sys.argv[4]))
