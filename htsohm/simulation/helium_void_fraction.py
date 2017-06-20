import sys
import os
import subprocess
import shutil
from datetime import datetime
from uuid import uuid4

import htsohm
from htsohm import config
from htsohm.material_files import write_cif_file, write_mixing_rules
from htsohm.material_files import write_pseudo_atoms, write_force_field
from htsohm.simulation.calculate_bin import calc_bin

def write_raspa_file(filename, uuid, simulation_config):
    """Writes RASPA input file for calculating helium void fraction.

    Args:
        filename (str): path to input file.
        run_id (str): identification string for run.
        material_id (str): uuid for material.

    Writes RASPA input-file.

    """
    simulation_cycles = simulation_config['simulation_cycles']
    with open(filename, "w") as raspa_input_file:
        raspa_input_file.write(
            "SimulationType         MonteCarlo\n" +
            "NumberOfCycles         %s\n" % simulation_cycles +     # number of MonteCarlo cycles
            "PrintEvery             10\n" +
            "PrintPropertiesEvery   10\n" +
            "\n" +
            "Forcefield             GenericMOFs\n" +
            "CutOff                 12.8\n" +           # LJ interaction cut-off, Angstroms
            "\n" +
            "Framework              0\n" +
            "FrameworkName          %s\n" % (uuid) +
            "UnitCells              1 1 1\n" +
            "ExternalTemperature    298.0\n" +    # External temperature, K
            "\n" +
            "Component 0 MoleculeName               helium\n" +
            "            MoleculeDefinition         TraPPE\n" +
            "            WidomProbability           1.0\n" +
            "            CreateNumberOfMolecules    0\n")

def parse_output(output_file, simulation_config):
    """Parse output file for void fraction data.

    Args:
        output_file (str): path to simulation output file.

    Returns:
        results (dict): average Widom Rosenbluth-weight.

    """
    results = {}
    with open(output_file) as origin:
        for line in origin:
            if not "Average Widom Rosenbluth-weight:" in line:
                continue
            results['vf_helium_void_fraction'] = float(line.split()[4])
        print("\nVOID FRACTION :   %s\n" % (results['vf_helium_void_fraction']))

    # calculate bin
    results['void_fraction_bin'] = calc_bin(
                results['vf_helium_void_fraction'],
                *simulation_config['limits'],
                config['number_of_convergence_bins'])

    return results

def run(material, simulation_config):
    """Runs void fraction simulation.

    Args:
        material (Material): material record.

    Returns:
        results (dict): void fraction simulation results.

    """
    simulation_directory  = config['simulations_directory']
    if simulation_directory == 'HTSOHM':
        htsohm_dir = os.path.dirname(os.path.dirname(htsohm.__file__))
        path = os.path.join(htsohm_dir, material.run_id)
    elif simulation_directory == 'SCRATCH':
        path = os.environ['SCRATCH']
    else:
        print('OUTPUT DIRECTORY NOT FOUND.')
    output_dir = os.path.join(path, 'output_%s_%s' % (material.uuid, uuid4()))
    print("Output directory :\t%s" % output_dir)
    os.makedirs(output_dir, exist_ok=True)
    filename = os.path.join(output_dir, "VoidFraction.input")
    write_raspa_file(filename, material.uuid, simulation_config)
    write_cif_file(material, output_dir)
    write_mixing_rules(material, output_dir)
    write_pseudo_atoms(material, output_dir)
    write_force_field(output_dir)
    while True:
        try:
            print("Date :\t%s" % datetime.now().date().isoformat())
            print("Time :\t%s" % datetime.now().time().isoformat())
            print("Calculating void fraction of %s..." % (material.uuid))
            subprocess.run(['simulate', './VoidFraction.input'], check=True, cwd=output_dir)
            filename = "output_%s_1.1.1_298.000000_0.data" % (material.uuid)
            output_file = os.path.join(output_dir, 'Output', 'System_0', filename)
            results = parse_output(output_file, simulation_config)
            shutil.rmtree(output_dir, ignore_errors=True)
            sys.stdout.flush()
        except (FileNotFoundError, IndexError, KeyError) as err:
            print(err)
            print(err.args)
            continue
        break

    return results
