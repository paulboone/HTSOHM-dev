
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

    bins = bng.check_number_of_bins(run_id)

    ml_step = ml_max / float(bins)
    sa_step = sa_max / float(bins)
    vf_step = vf_max / float(bins)

    ml_edges = np.arange( ml_min, ml_max + ml_step, ml_step )
    sa_edges = np.arange( sa_min, sa_max + sa_step, sa_step )
    vf_edges = np.arange( vf_min, vf_max + vf_step, vf_step )

    for i in range( bins ):
        if sa >= sa_edges[i] and sa < sa_edges[i + 1]:
            sa_bin = i
        if ml >= ml_edges[i] and ml < ml_edges[i + 1]:
            ml_bin = i
        if vf >= vf_edges[i] and vf < vf_edges[i + 1]:
            vf_bin = i

    run_data.methane_loading_bin = ml_bin
    run_data.surface_area_bin = sa_bin
    run_data.void_fraction_bin = vf_bin
    session.commit()

def run_all_simulations(id):
    run_data = session.query(RunData).get(id)
    
    ### RUN HELIUM VOID FRACTION
    results = helium_void_fraction_simulation.run(run_data.run_id, run_data.material_id)
    run_data.helium_void_fraction = results['VF_val']
    session.commit()
    
    ### RUN METHANE LOADING
    results = methane_loading_simulation.run(run_data.run_id, run_data.material_id, run_data.helium_void_fraction)
    run_data.absolute_volumetric_loading    = results['ML_a_cc']
    run_data.absolute_gravimetric_loading   = results['ML_a_cg']
    run_data.absolute_molar_loading         = results['ML_a_mk']
    run_data.excess_volumetric_loading      = results['ML_e_cc']
    run_data.excess_gravimetric_loading     = results['ML_e_cg']
    run_data.excess_molar_loading           = results['ML_e_mk']
    session.commit()
    
    ### RUN SURFACE AREA
    results = surface_area_simulation.run(run_data.run_id, run_data.material_id)
    run_data.unit_cell_surface_area     = results['SA_a2']
    run_data.volumetric_surface_area    = results['SA_mc']
    run_data.gravimetric_surface_area   = results['SA_mg']
    session.commit()
    
    # GetBins(id)

def dummy_test(run_id, generation):
#
#    import os
#    import subprocess
#    import shlex
#    import shutil
#
#    from sqlalchemy import create_engine
#    from sqlalchemy.orm import sessionmaker
#
#    from runDB_declarative import RunData, Base
#
#    import numpy as np

    tolerance = 0.05      # Acceptable deviation from original value(s)...
    number_of_trials = 1    # Number of times each simulation is repeated.

    wd = os.environ['HTSOHM_DIR']
    with open(wd + '/' + run_id + '.txt') as origin:
        for line in origin:
            if "Children per generation:" in line:
                children_per_generation = int(line.split()[3])
    first = generation * children_per_generation
    last = (generation + 1) * children_per_generation
    child_ids = np.arange(first, last)

    parents_ids = []
    for i in child_ids:

        parent_id = get_value(run_id, i, "parent_id")
        if parent_id not in parents_ids:
            parents_ids.append(parent_id)

    bad_materials = []
    for i in parents_ids:
        material_id = id_to_mat(run_id, i)

        print( "\nRe-Simulating %s-%s...\n" % (run_id, material_id) )

        ml_o = get_value(run_id, material_id, "absolute_volumetric_loading")
        sa_o = get_value(run_id, material_id, "volumetric_surface_area")
        vf_o = get_value(run_id, material_id, "helium_void_fraction")

        vfs = []
        for j in range(number_of_trials):
            void_fraction(run_id, material_id)
            vf_data = ( "Output/System_0/output_" +
                        "%s-%s_1.1.1_298.000000_0.data" % (run_id, material_id)
                        )
            with open(vf_data) as origin:
                for line in origin:
                    if not "Average Widom Rosenbluth-weight:" in line:
                        continue
                    try:
                        vf_val = line.split()[4]
                    except IndexError:
                        print()

            vfs.append( float(vf_val) )

        mls = []
        for j in range(number_of_trials):
            methane_loading(run_id, material_id)
            ml_data = ( "Output/System_0/output_" +
                        "%s-%s_1.1.1_298.000000_" % (run_id, material_id) +
                        "3.5e+06.data" )
            with open(ml_data) as origin:
                for line in origin:
                    if "absolute [cm^3 (STP)/c" in line:
                        ml_a_cc = line.split()[6]
            mls.append( float(ml_a_cc) )

        sas = []
        for j in range(number_of_trials):
            surface_area(run_id, material_id)
            sa_data = ( "Output/System_0/output_" +
                        "%s-%s_1.1.1_298.000000_0.data" % (run_id, material_id)
                        )
            with open(sa_data) as origin:
                count = 0
                for line in origin:
                    if "Surface area" in line:
                        if count == 0:
                            sa_a2 = line.split()[2]
                            count = count + 1
                        elif count == 1:
                            sa_mg = line.split()[2]
                            count = count + 1
                        elif count == 2:
                            sa_mc = line.split()[2]
            sas.append( float(sa_mc) )

        if abs(np.mean(mls) - ml_o) >= tolerance * ml_o:
            bad_materials.append(i)
        if abs(np.mean(sas) - sa_o) >= tolerance * sa_o:
            bad_materials.append(i)
        if abs(np.mean(vfs) - vf_o) >= tolerance * vf_o:
            bad_materials.append(i)

    failed = []
    for i in bad_materials:
        if i not in failed:
            failed.append(i)

    for i in parents_ids:
        material_id = id_to_mat(run_id, i)
        if i not in failed:
            data = {"dummy_test_result": 'y'}
            update_table(run_id, material_id, data)
        elif i in failed:
            data = {"dummy_test_result": 'n'}
            update_table(run_id, material_id, data)
#            AffectedMats = s.query(RunData).filter(RunData.run_id == run_id,
#                                                   RunData.parent_id == i
#                                                   ).all()
#            Maybes = [j.Mat for j in AffectedMats]
#            for j in Maybes:
#                update_table(run_id, j, data)

    if len(failed) == 0:
        print( "\nALL PARENTS IN GENERATION " +
               "%s PASSED THE DUMMY TEST.\n" % (generation) )
    if len(failed) != 0:
        print( "\nTHE FOLLOWING PARENTS IN GENERATION " +
               "%s FAIL THE DUMMY TEST:" % (generation) )
        for i in failed:
            material_id = id_to_mat(run_id, i)
            print( "\t%s-%s\n" % (run_id, material_id) )
