#! /usr/bin/env python

def simulate(run_ID, mat_ID):

    import os
    import subprocess
    import shlex
    import shutil
    
    from sqlalchemy import create_engine
    from sqlalchemy.orm import sessionmaker
    
    from runDB_declarative import RunData, Base

    import numpy as np

    print("\nMATERIAL ID :   %s-%s\n" % (run_ID, mat_ID) )
    pwd = os.getcwd()

    engine = create_engine( "sqlite:///HTSOHM-dev.db" )
    Base.metadata.bind = engine

    DBSession = sessionmaker(bind=engine)
    session = DBSession()

    check_first = session.query(RunData).filter( RunData.Run == run_ID,
                                                 RunData.Mat == str(mat_ID) ).count()
    if not check_first:
        new_mat = RunData( Run=run_ID, Mat=str(mat_ID) )
        session.add(new_mat)
        session.commit()

    DBmat = session.query(RunData).filter( RunData.Run == run_ID,
                                           RunData.Mat == str(mat_ID) )
#    new_mat = RunData(Run=run_ID, Mat=mat_ID)
#    session.add(new_mat)
#    session.commit()  

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
                    "            CreateNumberOfMolecules\t0" )
    VF_input.close()

    rd = os.environ['RASPA_DIR']

    subprocess.call(shlex.split( rd + '/bin/simulate -i VoidFraction.input' ))

    VF_data = "Output/System_0/output_%s-%s_1.1.1_298.000000_0.data" % (run_ID, mat_ID)
    with open(VF_data) as origin:
        for line in origin:
            if not "Average Widom Rosenbluth-weight:" in line:
                continue
            try:
                VF_val = line.split()[4]
            except IndexError:
                print()

    # Add to database...
    DBmat.update({'VF_wido': VF_val})
#    session.query(RunData).filter(RunData.Run == run_ID, RunData.Mat == mat_ID).update({'VF_wido': VF_val})
    session.commit()

    print( "\nVOID FRACTION :   %s\n" % (VF_val) )

    os.remove("VoidFraction.input")
    shutil.rmtree("Output")
    shutil.rmtree("Movies")
    shutil.rmtree("VTK")
    shutil.rmtree("Restart")

    # Simulate METHANE LOADING:
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
                    "HeliumVoidFraction %s\n" % (VF_val) +
                    "ExternalTemperature 298.0\n" +       # External temperature, K
                    "ExternalPressure 3500000\n" +        # External pressure, Pa
                    "\n" +
                    "Component 0 MoleculeName\t\tmethane\n" +
                    "            MoleculeDefinition\t\tTraPPE\n" +
                    "            TranslationProbability\t1.0\n" +
                    "            ReinsertionProbability\t1.0\n" +
                    "            SwapProbability\t\t1.0\n" +
                    "            CreateNumberOfMolecules\t0" )
    ML_input.close()

    subprocess.call(shlex.split( rd + '/bin/simulate -i MethaneLoading.input' ))

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

    # Add to database...
#    session.query(RunData).filter(RunData.Run == run_ID, RunData.Mat == mat_ID).update({'Abs_cccc': ML_a_cc, 'Abs_ccgr': ML_a_cg, 'Abs_mokg': ML_a_mk, 'Exc_cccc': ML_e_cc, 'Exc_ccgr': ML_e_cg, 'Exc_mokg': ML_e_mk})
    DBmat.update({'Abs_cccc': ML_a_cc, 'Abs_ccgr': ML_a_cg, 'Abs_mokg': ML_a_mk, 'Exc_cccc': ML_e_cc, 'Exc_ccgr': ML_e_cg, 'Exc_mokg': ML_e_mk})
    session.commit()

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
                    "            CreateNumberOfMolecules\t0" )
    SA_input.close()

    subprocess.call(shlex.split( rd + '/bin/simulate -i SurfaceArea.input' ))

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
    
    # Add to database...
#    session.query(RunData).filter(RunData.Run == run_ID, RunData.Mat == mat_ID).update({'SA_A2': SA_a2, 'SA_m2cc': SA_mc, 'SA_m2gr': SA_mg})
    DBmat.update({'SA_A2': SA_a2, 'SA_m2cc': SA_mc, 'SA_m2gr': SA_mg})
    session.commit()

    print( "\nSURFACE AREA\n" +
           "%s\tA^2\n" % (SA_a2) +
           "%s\tm^2/g\n" % (SA_mg) +
           "%s\tm^2/cm^3" % (SA_mc) )

    os.remove("SurfaceArea.input")
    shutil.rmtree("Output")
    shutil.rmtree("Movies")
    shutil.rmtree("VTK")
    shutil.rmtree("Restart")    

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

    SA = float(SA_mc)
    ML = float(ML_a_cc)
    VF = float(VF_val)

    for i in range( bins ):
        if SA >= SA_edges[i] and SA < SA_edges[i + 1]:
            SA_bin = i
        if ML >= ML_edges[i] and ML < ML_edges[i + 1]:
            ML_bin = i
        if VF >= VF_edges[i] and VF < VF_edges[i + 1]:
            VF_bin = i

    DBmat.update({'Bin_ML': str(ML_bin), 'Bin_SA': str(SA_bin), 'Bin_VF': str(VF_bin)})
    session.commit()

    
   
if __name__ == "__main__":
    import sys
    simulate(str(sys.argv[1]),
             int(sys.argv[2]))
