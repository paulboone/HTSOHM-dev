import os
import subprocess
import shutil

def write_raspa_file(filename, run_id, material_id):
    with open(filename, "w") as config:
        config.write("SimulationType\t\t\tMonteCarlo\n" +
                        "NumberOfCycles\t\t\t1000\n" +     # number of MonteCarlo cycles
                        "PrintEvery\t\t\t100\n" +
                        "PrintPropertiesEvery\t\t100\n" +
                        "\n" +
                        "Forcefield\t\t\t%s-%s\n" % (run_id, material_id) +
                        "CutOff\t\t\t\t12.8\n" +           # LJ interaction cut-off, Angstroms
                        "\n" +
                        "Framework 0\n" +
                        "FrameworkName %s-%s\n" % (run_id, material_id) +
                        "UnitCells 1 1 1\n" +
                        "ExternalTemperature 298.0\n" +    # External temperature, K
                        "\n" +
                        "Component 0 MoleculeName\t\thelium\n" +
                        "            MoleculeDefinition\t\tTraPPE\n" +
                        "            WidomProbability\t\t1.0\n" +
                        "            CreateNumberOfMolecules\t0\n" )

def parse_output(output_file):
    results = {}
    with open(output_file) as origin:
        for line in origin:
            if not "Average Widom Rosenbluth-weight:" in line:
                continue
            try:
                results['VF_val'] = float(line.split()[4][:-1])
            except IndexError:
                print()

    print( "\nVOID FRACTION :   %s\n" % (results['VF_val']) )

    return results

def run(run_id, material_id):
    os.makedirs('output', exist_ok=True)
    filename = "output/VoidFraction.input"
    write_raspa_file(filename, run_id, material_id)
    subprocess.run(['simulate', './VoidFraction.input'], check=True, cwd='output')

    output_file = "output/Output/System_0/output_%s-%s_1.1.1_298.000000_0.data" % (run_id, material_id)
    results = parse_output(output_file)
    shutil.rmtree("output")

    return results
