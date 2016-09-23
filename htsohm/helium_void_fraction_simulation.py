import os
import subprocess
import shutil

from htsohm import config

def write_raspa_file(filename, run_id, material_id):
    simulation_cycles = config['helium-void-fraction']['simulation-cycles']
    with open(filename, "w") as raspa_input_file:
        raspa_input_file.write(
            "SimulationType\t\t\tMonteCarlo\n" +
            "NumberOfCycles\t\t\t%s\n" % simulation_cycles +     # number of MonteCarlo cycles
            "PrintEvery\t\t\t10\n" +
            "PrintPropertiesEvery\t\t10\n" +
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
            "            CreateNumberOfMolecules\t0\n")

def parse_output(output_file):
    results = {}
    with open(output_file) as origin:
        for line in origin:
            if not "Average Widom Rosenbluth-weight:" in line:
                continue
            try:
                results['VF_val'] = float(line.split()[4])
            except IndexError:
                print()
    print("\nVOID FRACTION :   %s\n" % (results['VF_val']))

    return results

def run(run_id, material_id):
    simulation_directory  = config['simulations-directory']
    output_dir = os.path.join(os.environ[simulation_directory], 'output_%s' % material_id)
    os.makedirs(output_dir, exist_ok=True)
    filename = os.path.join(output_dir, "VoidFraction.input")
    write_raspa_file(filename, run_id, material_id)
    subprocess.run(['simulate', './VoidFraction.input'], check=True, cwd=output_dir)

    filename = "output_%s-%s_1.1.1_298.000000_0.data" % (run_id, material_id)
    output_file = os.path.join(output_dir, 'Output', 'System_0', filename)
    results = parse_output(output_file)
    shutil.rmtree(output_dir, ignore_errors=True)

    return results
