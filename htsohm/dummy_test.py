# standard library imports
#import os
#import shlex
#import shutil
#import subprocess

# related third party imports
import numpy as np

# local application/library specific imports
from htsohm import helium_void_fraction_simulation
from htsohm import methane_loading_simulation
from htsohm import surface_area_simulation
from htsohm.runDB_declarative import Material, session
#from htsohm.utilities import read_config_file

def dummy_test(run_id, id):
    """Recalculate material structure-properties to prevent statistical errors.

    Because methane loading, surface area, and helium void fractions are calculated using
    statistical methods (namely grand canonic Monte Carlo simulations) they are susceptible
    to statistical errors. To mitigate this, after a material has been selected as a potential
    parent, it's combination of structure-properties is resimulated some number of times and
    compared to the initally-calculated material. If the resimulated values differ from the
    initially-calculated value beyond an accpetable tolerance, the material fails the `dummy-test`
    and is flagged, preventing it from being used to generate new materials in the future.
    """
    tolerance = 0.5
    number_of_trials = 1

    material = session.query(Material).get(str(id))
    if material.dummy_test_result != 'pass':
        print( "\nRe-Simulating %s-%s...\n" % (run_id, id) )

        ####################################################################
        # resimulate helium void fraction(s)
        void_fractions = []
        for j in range(number_of_trials):
            vf_result = helium_void_fraction_simulation.run(run_id, id)
            void_fractions.append(float(vf_result["VF_val"]))

        ####################################################################
        # resimulate methane loading(s)
        methane_loadings = []
        for j in range(number_of_trials):
            ml_result = methane_loading_simulation.run(run_id, id, np.mean(void_fractions))
            methane_loadings.append(float(ml_result["ML_a_cc"]))

        ####################################################################
        # resimulate surface area(s)
        surface_areas = []
        for j in range(number_of_trials):
            sa_result = surface_area_simulation.run(run_id, id)
            surface_areas.append( float(sa_result["SA_mc"]) )

        ####################################################################
        # compare average of resimulated values to originals
        material = session.query(Material).get(id)
        ml_o = material.absolute_volumetric_loading    # initally-calculated values
        sa_o = material.volumetric_surface_area
        vf_o = material.helium_void_fraction
        if (
            abs(np.mean(methane_loadings) - ml_o) >= tolerance * ml_o or
            abs(np.mean(surface_areas) - sa_o) >= tolerance * sa_o or
            abs(np.mean(void_fractions) - vf_o) >= tolerance * vf_o
        ):
            material.dummy_test_result = 'fail'        # flag failed material
            print('%s-%s HAS FAILED PARENT-SCREENING.' % (run_id, id))
            return False        
        else:
            material.dummy_test_result = 'pass'        # prevents from future re-testing
            print('%s-%s HAS PASSED PARENT-SCREENING.' % (run_id, id))
            return True
    else:
        print('%s-%s HAD ALREADY PASSED PARENT-SCREENING.' % (run_id, id))
        return True
