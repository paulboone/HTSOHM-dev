import os
import subprocess
import shlex
import shutil

import numpy as np

from htsohm.runDB_declarative import Material, session
from htsohm import binning as bng
from htsohm import helium_void_fraction_simulation
from htsohm import methane_loading_simulation
from htsohm import surface_area_simulation
from htsohm.utilities import read_config_file

def get_bins(id, methane_loading, surface_area, void_fraction):
    """Returns methane_loading_bin, surface_area_bin, and void_fraction_bin.
    Each material is sorted into a bin corresponding to its combination of structure-properties.
    First, the structure property space is subdivided into arbitrary quadrants, or bins, then
    the simulated properties for a particular material are used to assigned it to a particular
    bin."""
    run_data = session.query(Material).get(id)

    ############################################################################
    # assign arbitrary maxima and subdivide the parameter space.
    config = read_config_file(run_data.run_id)
    bins = config["number-of-bins"]
    ml_min = 0.
    ml_max = 350.
    sa_min = 0.
    sa_max = 4500.
    vf_min = 0.
    vf_max = 1.
    ml_step = ml_max / float(bins)
    sa_step = sa_max / float(bins)
    vf_step = vf_max / float(bins)
    ml_edges = np.arange(ml_min, ml_max + ml_step, ml_step)
    sa_edges = np.arange(sa_min, sa_max + sa_step, sa_step)
    vf_edges = np.arange(vf_min, vf_max + vf_step, vf_step)

    ############################################################################
    # assign material to its respective bin
    for i in range( bins ):
        if surface_area >= sa_edges[i] and surface_area <= sa_edges[i + 1]:
            sa_bin = i
        if methane_loading >= ml_edges[i] and methane_loading <= ml_edges[i + 1]:
            ml_bin = i
        if void_fraction >= vf_edges[i] and void_fraction <= vf_edges[i + 1]:
            vf_bin = i
    print("\nBINS\t%s\t%s\t%s\n" % (ml_bin, sa_bin, vf_bin))

    results = {}
    results['ml_bin'] = ml_bin
    results['sa_bin'] = sa_bin
    results['vf_bin'] = vf_bin
    return results

def run_all_simulations(id):
    """Simulate helium void fraction, methane loading, and surface area.

    For a given material (id) three simulations are run using RASPA. First a helium void fraction
    is calculated, and then it is used to run a methane loading simulation (void fraction needed to
    calculate excess v. absolute loading). Finally, a surface area is calculated and the material is
    assigned to its appropriate bin."""
    run_data = session.query(Material).get(id)

    ############################################################################
    # run helium void fraction simulation
    results = helium_void_fraction_simulation.run(run_data.run_id, run_data.id)
    run_data.helium_void_fraction = results['VF_val']
    void_fraction = float(results['VF_val'])

    ############################################################################
    # run methane loading simulation
    results = methane_loading_simulation.run(run_data.run_id,
                                             run_data.id,
                                             run_data.helium_void_fraction)
    run_data.absolute_volumetric_loading   = results['ML_a_cc']
    run_data.absolute_gravimetric_loading  = results['ML_a_cg']
    run_data.absolute_molar_loading        = results['ML_a_mk']
    run_data.excess_volumetric_loading     = results['ML_e_cc']
    run_data.excess_gravimetric_loading    = results['ML_e_cg']
    run_data.excess_molar_loading          = results['ML_e_mk']
    run_data.host_host_avg                 = results['host_host_avg']
    run_data.host_host_vdw                 = results['host_host_vdw']
    run_data.host_host_cou                 = results['host_host_cou']
    run_data.adsorbate_adsorbate_avg       = results['adsorbate_adsorbate_avg']
    run_data.adsorbate_adsorbate_vdw       = results['adsorbate_adsorbate_vdw']
    run_data.adsorbate_adsorbate_cou       = results['adsorbate_adsorbate_cou']
    run_data.host_adsorbate_avg            = results['host_adsorbate_avg']
    run_data.host_adsorbate_vdw            = results['host_adsorbate_vdw']
    run_data.host_adsorbate_cou            = results['host_adsorbate_cou']
    methane_loading = float(results['ML_a_cc'])

    ############################################################################
    # run surface area simulation
    results = surface_area_simulation.run(run_data.run_id, run_data.id)
    run_data.unit_cell_surface_area     = results['SA_a2']
    run_data.volumetric_surface_area    = results['SA_mc']
    run_data.gravimetric_surface_area   = results['SA_mg']
    surface_area = float(results['SA_mc'])

    ############################################################################
    # assign material to bin
    results = get_bins(run_data.id, methane_loading, surface_area, void_fraction)
    run_data.methane_loading_bin = results['ml_bin']
    run_data.surface_area_bin = results['sa_bin']
    run_data.void_fraction_bin = results['vf_bin']

    run_data.data_complete = True

def dummy_test(run_id, next_generation, generation):
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

    failed = []
    for material in next_generation:
        ########################################################################
        # iterate over all selected-parents for the next generation
        parent = session.query(Material).get(material.parent_id)

        if parent.dummy_test_result != "pass":       # materials are not re-tested
            print( "\nRe-Simulating %s-%s...\n" % (run_id, parent.id) )

            ####################################################################
            # resimulate helium void fraction(s)
            void_fractions = []
            for j in range(number_of_trials):
                vf_result = helium_void_fraction_simulation.run(run_id, parent.id)
                void_fractions.append(float(vf_result["VF_val"]))

            ####################################################################
            # resimulate methane loading(s)
            methane_loadings = []
            for j in range(number_of_trials):
                ml_result = methane_loading_simulation.run(run_id, parent.id,
                    np.mean(void_fractions))
                methane_loadings.append(float(ml_result["ML_a_cc"]))

            ####################################################################
            # resimulate surface area(s)
            surface_areas = []
            for j in range(number_of_trials):
                sa_result = surface_area_simulation.run(run_id, parent.id)
                surface_areas.append( float(sa_result["SA_mc"]) )

            ####################################################################
            # compare average of resimulated values to originals
            ml_o = parent.absolute_volumetric_loading    # initally-calculated values
            sa_o = parent.volumetric_surface_area
            vf_o = parent.helium_void_fraction
            if ( abs(np.mean(methane_loadings) - ml_o) >= tolerance * ml_o or
                 abs(np.mean(surface_areas) - sa_o) >= tolerance * sa_o or
                 abs(np.mean(void_fractions) - vf_o) >= tolerance * vf_o ):
                parent.dummy_test_result = "fail"        # flag failed material
                print(
                    "A MATERIAL HAS FAILED!\n" +
                    "Run:\t%s\n" % (run_id) +
                    "Material:\t%s\n" % (parent.id))
                failed.append("%s-%s" % (run_id, parent.id))
                break
            else:
                parent.dummy_test_result = "pass"        # prevents from future re-testing
    ############################################################################
    # if any materials fail, then new parents are selected and tested (see : htsohm/htsohm.py)
    if len(failed) == 0:
        return True # Success!

    return False #Fail
