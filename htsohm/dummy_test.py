# standard library imports
#import os
#import shlex
#import shutil
#import subprocess

# related third party imports
import numpy as np

# local application/library specific imports
from htsohm.binning import select_parent
from htsohm.simulate import run_all_simulations
from htsohm.utilities import read_config_file
from htsohm import helium_void_fraction_simulation
from htsohm import methane_loading_simulation
from htsohm import surface_area_simulation
from htsohm.db import Material, session

def screen_parent(run_id):
    test_complete = False
    while not test_complete:
        parent_id = select_parent(run_id)
        test_complete = dummy_test(run_id, parent_id)
    return parent_id

def retest(m_orig, retests, tolerance):
    """Recalculate material structure-properties to prevent statistical errors.

    Because methane loading, surface area, and helium void fractions are calculated using
    statistical methods (namely grand canonic Monte Carlo simulations) they are susceptible
    to statistical errors. To mitigate this, after a material has been selected as a potential
    parent, it's combination of structure-properties is resimulated some number of times and
    compared to the initally-calculated material. If the resimulated values differ from the
    initially-calculated value beyond an accpetable tolerance, the material fails the `dummy-test`
    and is flagged, preventing it from being used to generate new materials in the future.
    """

    m = m_orig.clone()
    run_all_simulations(m)

    # requery row from database, in case someone else has changed it, and lock it
    # if the row is presently locked, this method blocks until the row lock is released
    session.refresh(m_orig, lockmode='update')
    if m_orig.retest_num < retests:
        m_orig.retest_methane_loading_sum += m.absolute_volumetric_loading
        m_orig.retest_surface_area_sum += m.volumetric_surface_area
        m_orig.retest_void_fraction_sum += m.helium_void_fraction
        m_orig.retest_num += 1

        if m_orig.retest_num == retests:
            m_orig.retest_passed = m.calculate_retest_result(tolerance)
    else:
        pass
        # otherwise our test is extra / redundant and we don't save it

    session.commit()
