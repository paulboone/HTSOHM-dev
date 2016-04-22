#! /usr/bin/env python

import os
import subprocess
import shlex
import shutil

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
import numpy as np

from htsohm.runDB_declarative import Base, RunData, CreateSession

def CreateSession():
#
#    from sqlalchemy import create_engine
#    from sqlalchemy.orm import sessionmaker
#    from runDB_declarative import Base

    engine = create_engine( "sqlite:///HTSOHM-dev.db" )
    Base.metadata.bind = engine

    DBSession = sessionmaker(bind=engine)
    session = DBSession()

    return session

def AddRows(run_ID, mat_IDs):
#
#    from runDB_declarative import RunData

    s = CreateSession()

    for i in mat_IDs:
        check_first = s.query(RunData).filter( RunData.run_id == run_ID,
                                               RunData.material_id == str(i)
                                              ).count()
        if not check_first:
            new_mat = RunData( run_id=run_ID, material_id=str(i) )
            s.add(new_mat)
    s.commit()


def UpdateTable(run_ID, mat_ID, data):
#
#   from runDB_declarative import RunData

    s = CreateSession()
    DBmat = s.query(RunData).filter( RunData.run_id == run_ID,
                                     RunData.material_id == str(mat_ID) )
    DBmat.update(data)
    s.commit()


def GetValue(run_ID, mat_ID, value):
#
#    from runDB_declarative import RunData

    s = CreateSession()
    DBmat = s.query(RunData).filter( RunData.run_id == run_ID,
                                     RunData.material_id == str(mat_ID) )

    for i in DBmat:
        DBval = getattr(i, value)

    return DBval

def id_to_mat(run_ID, ID):
#
#    from runDB_declarative import RunData

    s = CreateSession()
    DBmat = s.query(RunData).filter( RunData.run_id == run_ID,
                                     RunData.id == str(ID) )

    for i in DBmat:
        mat = getattr(i, "material_id")

    return mat


def VoidFraction(run_ID, mat_ID):
#
#    import os
#    import subprocess
#    import shlex

    pwd = os.getcwd()
    
    # Simulate VOID FRACTION:
    VF_input = open( pwd + '/VoidFraction.input', "w")
    VF_input.write( "SimulationType\t\t\tMonteCarlo\n" +
                    "NumberOfCycles\t\t\t1000\n" +             # number of MonteCarlo cycles
         "PrintEvery\t\t\t100\n" +
                    "PrintPropertiesEvery\t\t100\n" +
                    "\n" +
                    "Forcefield\t\t\t%s-%s\n" % (run_ID, mat_ID) +
                    "CutOff\t\t\t\t12.8\n" +                       # LJ interaction cut-off, Angstroms
                    "\n" +
                    "Framework 0\n" +
                    "FrameworkName %s-%s\n" % (run_ID, mat_ID) +
                    "UnitCells 1 1 1\n" +
                    "ExternalTemperature 298.0\n" +       # External temperature, K
                    "\n" +
                    "Component 0 MoleculeName\t\thelium\n" +
                    "            MoleculeDefinition\t\tTraPPE\n" +
                    "            WidomProbability\t\t1.0\n" +
                    "            CreateNumberOfMolecules\t0\n" )
    VF_input.close()

    subprocess.call(shlex.split('simulate VoidFraction.input'))

def MethaneLoading(run_ID, mat_ID):
#
#    import os
#    import subprocess
#    import shlex

    pwd = os.getcwd()

    VF = GetValue(run_ID, mat_ID, "helium_void_fraction")

    # Simulate METHANE LOADING
    ML_input = open( pwd + '/MethaneLoading.input', "w")
    ML_input.write( "SimulationType\t\t\tMonteCarlo\n" +
                    "NumberOfCycles\t\t\t1000\n" +             # number of MonteCarlo cycles
                    "NumberOfInitializationCycles\t500\n" +    # number of initialization cycles
                    "PrintEvery\t\t\t100\n" +
                    "RestartFile\t\t\tno\n" +
                    "\n" +
                    "Forcefield\t\t\t%s-%s\n" % (run_ID, mat_ID) +
                    "ChargeMethod\t\t\tEwald\n"
                    "CutOff\t\t\t\t12.0\n" +                   # electrostatic cut-off, Angstroms
                    "\n" +
                    "Framework 0\n" +
                    "FrameworkName %s-%s\n" % (run_ID, mat_ID) +
                    "UnitCells 1 1 1\n" +
                    "HeliumVoidFraction %s\n" % (VF) +
                    "ExternalTemperature 298.0\n" +            # External temperature, K
                    "ExternalPressure 3500000\n" +             # External pressure, Pa
                    "\n" +
                    "Component 0 MoleculeName\t\tmethane\n" +
                    "            MoleculeDefinition\t\tTraPPE\n" +
                    "            TranslationProbability\t1.0\n" +
                    "            ReinsertionProbability\t1.0\n" +
                    "            SwapProbability\t\t1.0\n" +
                    "            CreateNumberOfMolecules\t0\n" )
    ML_input.close()

    subprocess.call(shlex.split('simulate MethaneLoading.input'))


def SurfaceArea(run_ID, mat_ID):
#
#    import os
#    import subprocess
#    import shlex

    pwd = os.getcwd()

    # Simulate SURFACE AREA:
    SA_input = open( pwd + '/SurfaceArea.input', "w")
    SA_input.write( "SimulationType\t\t\tMonteCarlo\n" +
                    "NumberOfCycles\t\t\t10\n" +             # number of MonteCarlo cycles
                    "PrintEvery\t\t\t1\n" +
                    "PrintPropertiesEvery\t\t1\n" +
                    "\n" +
                    "Forcefield %s-%s\n" % (run_ID, mat_ID) +
                    "CutOff 12.8\n" +                        # electrostatic cut-off, Angstroms
                    "\n" +
                    "Framework 0\n" +
                    "FrameworkName %s-%s\n" % (run_ID, mat_ID) +
                    "UnitCells 1 1 1\n" +
                    "SurfaceAreaProbeDistance Minimum\n" +
                    "\n" +
                    "Component 0 MoleculeName\t\tN2\n" +
                    "            StartingBead\t\t0\n" +
                    "            MoleculeDefinition\t\tTraPPE\n" +
                    "            SurfaceAreaProbability\t1.0\n" +
                    "            CreateNumberOfMolecules\t0\n" )
    SA_input.close()

    subprocess.call(shlex.split('simulate SurfaceArea.input'))


def GetML(run_ID, mat_ID):

    MethaneLoading(run_ID, mat_ID)

    ML_data = "Output/System_0/output_%s-%s_1.1.1_298.000000_3.5e+06.data" % (run_ID, mat_ID)
    with open(ML_data) as origin:
        for line in origin:
            if "absolute [mol/kg" in line:
                ML_a_mk = line.split()[5]
            elif "absolute [cm^3 (STP)/g" in line:
                ML_a_cg = line.split()[6]
            elif "absolute [cm^3 (STP)/c" in line:
                ML_a_cc = line.split()[6]
            elif "excess [mol/kg" in line:
                ML_e_mk = line.split()[5]
            elif "excess [cm^3 (STP)/g" in line:
                ML_e_cg = line.split()[6]
            elif "excess [cm^3 (STP)/c" in line:
                ML_e_cc = line.split()[6]

    data = {"absolute_volumetric_loading": ML_a_cc,
            "absolute_gravimetric_loading": ML_a_cg,
            "absolute_molar_loading": ML_a_mk,
            "excess_volumetric_loading": ML_e_cc,
            "excess_gravimetric_loading": ML_e_cg,
            "excess_molar_loading": ML_e_mk}
    UpdateTable(run_ID, mat_ID, data)

    print( "\nMETHANE LOADING\tabsolute\texcess\n" +
           "mol/kg\t\t%s\t%s\n" % (ML_a_mk, ML_e_mk) +
           "cc/g\t\t%s\t%s\n" % (ML_a_cg, ML_e_cg) +
           "cc/cc\t\t%s\t%s\n" % (ML_a_cc, ML_e_cc) )

    #STILL NEED TO GREP HEATDESORP
    os.remove("MethaneLoading.input")
    shutil.rmtree("Output")
    shutil.rmtree("Movies")
    shutil.rmtree("VTK")
    shutil.rmtree("Restart")


def GetSA(run_ID, mat_ID):
    SurfaceArea(run_ID, mat_ID)

    SA_data = "Output/System_0/output_%s-%s_1.1.1_298.000000_0.data" % (run_ID, mat_ID)
    with open(SA_data) as origin:
        count = 0
        for line in origin:
            if "Surface area" in line:
                if count == 0:
                    SA_a2 = line.split()[2]
                    count = count + 1
                elif count == 1:
                    SA_mg = line.split()[2]
                    count = count + 1
                elif count == 2:
                    SA_mc = line.split()[2]
    
    data = {"unit_cell_surface_area": SA_a2,
            "volumetric_surface_area": SA_mc,
            "gravimetric_surface_area": SA_mg}
    UpdateTable(run_ID, mat_ID, data)

    print( "\nSURFACE AREA\n" +
           "%s\tA^2\n" % (SA_a2) +
           "%s\tm^2/g\n" % (SA_mg) +
           "%s\tm^2/cm^3" % (SA_mc) )

    os.remove("SurfaceArea.input")
    shutil.rmtree("Output")
    shutil.rmtree("Movies")
    shutil.rmtree("VTK")
    shutil.rmtree("Restart")


def GetVF(run_ID, mat_ID):

    VoidFraction(run_ID, mat_ID)

    VF_data = "Output/System_0/output_%s-%s_1.1.1_298.000000_0.data" % (run_ID, mat_ID)
    with open(VF_data) as origin:
        for line in origin:
            if not "Average Widom Rosenbluth-weight:" in line:
                continue
            try:
                VF_val = float(line.split()[4][:-1])
            except IndexError:
                print()

    print( "\nVOID FRACTION :   %s\n" % (VF_val) )

    # Add to database...
    data = {"helium_void_fraction": VF_val}
    UpdateTable(run_ID, mat_ID, data)

    os.remove("VoidFraction.input")
    shutil.rmtree("Output")
    shutil.rmtree("Movies")
    shutil.rmtree("VTK")
    shutil.rmtree("Restart")


def GetBins(run_ID, mat_ID):
   
    ML = GetValue(run_ID, mat_ID, "absolute_volumetric_loading")
    SA = GetValue(run_ID, mat_ID, "volumetric_surface_area")
    VF = GetValue(run_ID, mat_ID, "helium_void_fraction")
 
    # Arbitary structure-property space "boundaries"
    ML_min = 0.
    ML_max = 350.
    SA_min = 0.
    SA_max = 4500.
    VF_min = 0.
    VF_max = 1.

    bins = int( run_ID[-1] )

    ML_step = ML_max / float(bins)
    SA_step = SA_max / float(bins)
    VF_step = VF_max / float(bins)

    ML_edges = np.arange( ML_min, ML_max + ML_step, ML_step )
    SA_edges = np.arange( SA_min, SA_max + SA_step, SA_step )
    VF_edges = np.arange( VF_min, VF_max + VF_step, VF_step )

    for i in range( bins ):
        if SA >= SA_edges[i] and SA < SA_edges[i + 1]:
            SA_bin = i
        if ML >= ML_edges[i] and ML < ML_edges[i + 1]:
            ML_bin = i
        if VF >= VF_edges[i] and VF < VF_edges[i + 1]:
            VF_bin = i

    data = {"methane_loading_bin": str(ML_bin),
            "surface_area_bin": str(SA_bin),
            "void_fraction_bin": str(VF_bin)}
    UpdateTable(run_ID, mat_ID, data)


def simulate(run_ID, mat_IDs):
   for i in mat_IDs:
       GetVF(run_ID, i)
       GetML(run_ID, i)
       GetSA(run_ID, i)
       GetBins(run_ID, i)


def DummyTest(run_ID, generation):
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

    s = CreateSession()

    Tolerance = 0.05      # Acceptable deviation from original value(s)...
    NumberOfTrials = 1    # Number of times each simulation is repeated.

    wd = os.environ['HTSOHM_DIR']
    with open(wd + '/' + run_ID + '.txt') as origin:
        for line in origin:
            if "Children per generation:" in line:
                children_per_generation = int(line.split()[3])
    First = generation * children_per_generation
    Last = (generation + 1) * children_per_generation
    c_IDs = np.arange(First, Last)

    p_IDs = []
    for i in c_IDs:

        pID = GetValue(run_ID, i, "parent_id")
        if pID not in p_IDs:
            p_IDs.append(pID)

    BadMaterials = []
    for i in p_IDs:
        matID = id_to_mat(run_ID, i)

        print( "\nRe-Simulating %s-%s...\n" % (run_ID, matID) )

        ML_o = GetValue(run_ID, matID, "absolute_volumetric_loading")
        SA_o = GetValue(run_ID, matID, "volumetric_surface_area")
        VF_o = GetValue(run_ID, matID, "helium_void_fraction")

        VFs = []
        for j in range(NumberOfTrials):
            VoidFraction(run_ID, matID)
            VF_data = "Output/System_0/output_%s-%s_1.1.1_298.000000_0.data" % (run_ID, matID)
            with open(VF_data) as origin:
                for line in origin:
#                    if not "Average Widom:" in line:
#                        continue
#                    try:
#                        VF_val = line.split()[3]
#                    except IndexError:
#                        print()
                    if not "Average Widom Rosenbluth-weight:" in line:       #ONLY ON AKAIJA BUILD
                        continue
                    try:
                        VF_val = line.split()[4]
                    except IndexError:
                        print()

            VFs.append( float(VF_val) )

        MLs = []
        for j in range(NumberOfTrials):
            MethaneLoading(run_ID, matID)
            ML_data = "Output/System_0/output_%s-%s_1.1.1_298.000000_3.5e+06.data" % (run_ID, matID)
            with open(ML_data) as origin:
                for line in origin:
                    if "absolute [cm^3 (STP)/c" in line:
                        ML_a_cc = line.split()[6]
            MLs.append( float(ML_a_cc) )

        SAs = []
        for j in range(NumberOfTrials):
            SurfaceArea(run_ID, matID)
            SA_data = "Output/System_0/output_%s-%s_1.1.1_298.000000_0.data" % (run_ID, matID)
            with open(SA_data) as origin:
                count = 0
                for line in origin:
                    if "Surface area" in line:
                        if count == 0:
                            SA_a2 = line.split()[2]
                            count = count + 1
                        elif count == 1:
                            SA_mg = line.split()[2]
                            count = count + 1
                        elif count == 2:
                            SA_mc = line.split()[2]
            SAs.append( float(SA_mc) )

        if abs(np.mean(MLs) - ML_o) >= Tolerance * ML_o:
            BadMaterials.append(i)
        if abs(np.mean(SAs) - SA_o) >= Tolerance * SA_o:
            BadMaterials.append(i)
        if abs(np.mean(VFs) - VF_o) >= Tolerance * VF_o:
            BadMaterials.append(i)

    Failed = []
    for i in BadMaterials:
        if i not in Failed:
            Failed.append(i)

    for i in p_IDs:
        matID = id_to_mat(run_ID, i)
        if i not in Failed:
            data = {"dummy_test_result": 'y'}
            UpdateTable(run_ID, matID, data)
        elif i in Failed:
            data = {"dummy_test_result": 'n'}
            UpdateTable(run_ID, matID, data)
#            AffectedMats = s.query(RunData).filter(RunData.run_id == run_ID,
#                                                   RunData.parent_id == i
#                                                   ).all()
#            Maybes = [j.Mat for j in AffectedMats]
#            for j in Maybes:
#                UpdateTable(run_ID, j, data)
                                                    
    if len(Failed) == 0:
        print( "\nALL PARENTS IN GENERATION %s PASSED THE DUMMY TEST.\n" % (generation) )
    if len(Failed) != 0:
        print( "\nTHE FOLLOWING PARENTS IN GENERATION %s FAIL THE DUMMY TEST:" % (generation) )
        for i in Failed:
            matID = id_to_mat(run_ID, i)
            print( "\t%s-%s\n" % (run_ID, matID) )


if __name__ == "__main__":
    import sys
    simulate(str(sys.argv[1]),
             int(sys.argv[2]))
