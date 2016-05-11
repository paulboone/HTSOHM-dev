import os
import subprocess
import shutil

def write_raspa_file(filename, run_id, material_id):
    with open(filename, "w") as config:
        config.write(
            "SimulationType\t\t\tMonteCarlo\n" +
            "NumberOfCycles\t\t\t10\n" +             # number of MonteCarlo cycles
            "PrintEvery\t\t\t1\n" +
            "PrintPropertiesEvery\t\t1\n" +
            "\n" +
            "Forcefield %s-%s\n" % (run_id, material_id) +
            "CutOff 12.8\n" +                        # electrostatic cut-off, Angstroms
            "\n" +
            "Framework 0\n" +
            "FrameworkName %s-%s\n" % (run_id, material_id) +
            "UnitCells 1 1 1\n" +
            "SurfaceAreaProbeDistance Minimum\n" +
            "\n" +
            "Component 0 MoleculeName\t\tN2\n" +
            "            StartingBead\t\t0\n" +
            "            MoleculeDefinition\t\tTraPPE\n" +
            "            SurfaceAreaProbability\t1.0\n" +
            "            CreateNumberOfMolecules\t0\n")

def parse_output(output_file):
    results = {}
    with open(output_file) as origin:
        count = 0
        for line in origin:
            if "Surface area" in line:
                if count == 0:
                    results['SA_a2'] = line.split()[2]
                    count = count + 1
                elif count == 1:
                    results['SA_mg'] = line.split()[2]
                    count = count + 1
                elif count == 2:
                    results['SA_mc'] = line.split()[2]
    
    print(
        "\nSURFACE AREA\n" +
        "%s\tA^2\n"      % (results['SA_a2']) +
        "%s\tm^2/g\n"    % (results['SA_mg']) +
        "%s\tm^2/cm^3"   % (results['SA_mc']))
    
    return results

def run(run_id, material_id):
    os.makedirs('output', exist_ok=True)
    filename = 'output/SurfaceArea.input'
    write_raspa_file(filename, run_id, material_id)
    subprocess.run(
        ['simulate', './SurfaceArea.input'], check=True, cwd='output')

    output_dir = "output/Output/System_0/"
    output_file = "%soutput_%s-%s_1.1.1_298.000000_0.data" % (
        output_dir, run_id, material_id)
    results = parse_output(output_file)
    shutil.rmtree("output")

    return results
