import os
import subprocess
import shutil

def write_raspa_file(filename, run_id, material_id, helium_void_fraction ):
    with open(filename, "w") as config:
        config.write("SimulationType\t\t\tMonteCarlo\n" +
                        "NumberOfCycles\t\t\t1000\n" +             # number of MonteCarlo cycles
                        "NumberOfInitializationCycles\t500\n" +    # number of initialization cycles
                        "PrintEvery\t\t\t100\n" +
                        "RestartFile\t\t\tno\n" +
                        "\n" +
                        "Forcefield\t\t\t%s-%s\n" % (run_id, material_id) +
                        "ChargeMethod\t\t\tEwald\n"
                        "CutOff\t\t\t\t12.0\n" +                   # electrostatic cut-off, Angstroms
                        "\n" +
                        "Framework 0\n" +
                        "FrameworkName %s-%s\n" % (run_id, material_id) +
                        "UnitCells 1 1 1\n" +
                        "HeliumVoidFraction %s\n" % (helium_void_fraction) +
                        "ExternalTemperature 298.0\n" +            # External temperature, K
                        "ExternalPressure 3500000\n" +             # External pressure, Pa
                        "\n" +
                        "Component 0 MoleculeName\t\tmethane\n" +
                        "            MoleculeDefinition\t\tTraPPE\n" +
                        "            TranslationProbability\t1.0\n" +
                        "            ReinsertionProbability\t1.0\n" +
                        "            SwapProbability\t\t1.0\n" +
                        "            CreateNumberOfMolecules\t0\n" )

def parse_output(output_file):
    results = {}
    with open(output_file) as origin:
        for line in origin:
            if "absolute [mol/kg" in line:
                results['ML_a_mk'] = line.split()[5]
            elif "absolute [cm^3 (STP)/g" in line:
                results['ML_a_cg'] = line.split()[6]
            elif "absolute [cm^3 (STP)/c" in line:
                results['ML_a_cc'] = line.split()[6]
            elif "excess [mol/kg" in line:
                results['ML_e_mk'] = line.split()[5]
            elif "excess [cm^3 (STP)/g" in line:
                results['ML_e_cg'] = line.split()[6]
            elif "excess [cm^3 (STP)/c" in line:
                results['ML_e_cc'] = line.split()[6]
    
    print("\nMETHANE LOADING\tabsolute\texcess\n" +
            "mol/kg\t\t%s\t%s\n" % (results['ML_a_mk'], results['ML_e_mk']) +
            "cc/g\t\t%s\t%s\n"   % (results['ML_a_cg'], results['ML_e_cg']) +
            "cc/cc\t\t%s\t%s\n"  % (results['ML_a_cc'], results['ML_e_cc']) )
    
    return results

def run(run_id, material_id, helium_void_fraction):
    # run_data = session.query(RunData).get(id)
    
    os.makedirs('output', exist_ok=True)
    filename = 'output/MethaneLoading.input'
    write_raspa_file(filename, run_id, material_id, helium_void_fraction)
    subprocess.run(['simulate', './MethaneLoading.input'], check=True, cwd='output')
    
    output_file = "output/Output/System_0/output_%s-%s_1.1.1_298.000000_3.5e+06.data" % (run_id, material_id)
    results = parse_output(output_file)

    #STILL NEED TO GREP HEATDESORP
    shutil.rmtree("output")
    
    return results
