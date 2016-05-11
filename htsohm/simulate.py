import os
import subprocess
import shlex
import shutil

import numpy as np

from htsohm.runDB_declarative import RunData, session
from htsohm import binning as bng
from htsohm import helium_void_fraction_simulation
from htsohm import methane_loading_simulation
from htsohm import surface_area_simulation

def id_to_mat(id):
    return session.query(RunData).get(id).material_id

def get_bins(id):
    run_data = session.query(RunData).get(id)
    ml = run_data.absolute_volumetric_loading
    sa = run_data.volumetric_surface_area
    vf = run_data.helium_void_fraction

    # Arbitary structure-property space "boundaries"
    ml_min = 0.
    ml_max = 350.
    sa_min = 0.
    sa_max = 4500.
    vf_min = 0.
    vf_max = 1.

    bins = bng.check_number_of_bins(run_data.run_id)

    ml_step = ml_max / float(bins)
    sa_step = sa_max / float(bins)
    vf_step = vf_max / float(bins)

    ml_edges = np.arange(ml_min, ml_max + ml_step, ml_step)
    sa_edges = np.arange(sa_min, sa_max + sa_step, sa_step)
    vf_edges = np.arange(vf_min, vf_max + vf_step, vf_step)

    for i in range( bins ):
        if sa >= sa_edges[i] and sa < sa_edges[i + 1]:
            sa_bin = i
        if ml >= ml_edges[i] and ml < ml_edges[i + 1]:
            ml_bin = i
        if vf >= vf_edges[i] and vf < vf_edges[i + 1]:
            vf_bin = i

    print("\nBINS\t%s\t%s\t%s\n" % (ml_bin, sa_bin, vf_bin))
    results = {}
    results['ml_bin'] = ml_bin
    results['sa_bin'] = sa_bin
    results['vf_bin'] = vf_bin

    return results

def run_all_simulations(id):
    run_data = session.query(RunData).get(id)
    
    ### RUN HELIUM VOID FRACTION
    results = helium_void_fraction_simulation.run(run_data.run_id, run_data.material_id)
    run_data.helium_void_fraction = results['VF_val']
    session.commit()
    
    ### RUN METHANE LOADING
    results = methane_loading_simulation.run(run_data.run_id, 
                                             run_data.material_id,
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
    session.commit()
    
    ### RUN SURFACE AREA
    results = surface_area_simulation.run(run_data.run_id, run_data.material_id)
    run_data.unit_cell_surface_area     = results['SA_a2']
    run_data.volumetric_surface_area    = results['SA_mc']
    run_data.gravimetric_surface_area   = results['SA_mg']
    session.commit()
    
    ### GET BIN
    results = get_bins(id)
    run_data.methane_loading_bin = results['ml_bin']
    run_data.surface_area_bin = results['sa_bin']
    run_data.void_fraction_bin = results['vf_bin']
    session.commit()

def dummy_test(run_id, next_materials_list, status, generation):
    tolerance = 0.5
    number_of_trials = 1
    failed = []
    for i in next_materials_list:
        parent_id = str(i[1])
        parent = session.query(RunData).get(parent_id)
        
        if parent.dummy_test_result != "pass":
            print( "\nRe-Simulating %s-%s...\n" % (run_id, parent.material_id) )

            void_fractions = []                   # re-simulate void fraction calculations
            for j in range(number_of_trials):
                vf_result = helium_void_fraction_simulation.run(
                    run_id, parent.material_id)
                void_fractions.append( float(vf_result["VF_val"]) )

            methane_loadings = []                   # re-simulate methane loading calculations
            for j in range(number_of_trials):
                ml_result = methane_loading_simulation.run(
                    run_id, parent.material_id, np.mean(void_fractions))
                methane_loadings.append( float(ml_result["ML_a_cc"]) )

            surface_areas = []                   # re-simulate surface area calculations
            for j in range(number_of_trials):
                sa_result = surface_area_simulation.run(
                    run_id, parent.material_id)
                surface_areas.append( float(sa_result["SA_mc"]) )

            ml_o = parent.absolute_volumetric_loading
            sa_o = parent.volumetric_surface_area
            vf_o = parent.helium_void_fraction
            
            if ( abs(np.mean(methane_loadings) - ml_o) >= tolerance * ml_o or
                 abs(np.mean(surface_areas) - sa_o) >= tolerance * sa_o or
                 abs(np.mean(void_fractions) - vf_o) >= tolerance * vf_o ):
                parent.dummy_test_result = "fail"
                print( 
                    "A MATERIAL HAS FAILED!\n" +
                    "Run:\t%s\n" % (run_id) +
                    "Material:\t%s\n" % (parent.material_id))
                failed.append("%s-%s" % (run_id, parent.material_id))
                break
            else:
                parent.dummy_test_result = "pass"
            session.commit()

    if len(failed) == 0:
        status = "Dummy test:   COMPLETE"
    print (status)

    return status
