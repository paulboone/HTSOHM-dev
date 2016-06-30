# standard library imports
import os
import shlex
import shutil
import subprocess

# related third party imports
import numpy as np

# local application/library specific imports
from htsohm import helium_void_fraction_simulation
from htsohm import methane_loading_simulation
from htsohm import surface_area_simulation
from htsohm.runDB_declarative import Material, session
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
