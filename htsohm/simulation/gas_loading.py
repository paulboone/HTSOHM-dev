from datetime import datetime
from glob import glob
import math
import os
from pathlib import Path
import subprocess
import shutil
from string import Template
import sys
from uuid import uuid4

import numpy as np

from htsohm.material_files import write_mol_file, write_mixing_rules
from htsohm.material_files import write_pseudo_atoms, write_force_field
from htsohm.simulation.templates import load_and_subs_template
from htsohm.db import GasLoading

def write_raspa_file(filename, material, simulation_config, restart):
    """Writes RASPA input file for simulating gas adsorption.

    Args:
        filename (str): path to input file.
        material (Material): material record.

    Writes RASPA input-file.

    """
    # Set default void fraction value if non was calculated
    if len(material.void_fraction) == 0 or material.void_fraction[0].void_fraction is None:
        void_fraction = 0.
    else:
        void_fraction = material.void_fraction[0].void_fraction

    # Load simulation parameters from config
    unit_cells = material.structure.minimum_unit_cells(simulation_config['cutoff'])
    values = {
            "Restart"                       : 'yes' if restart else 'no',
            "Cutoff"                        : simulation_config['cutoff'],
            "NumberOfCycles"                : simulation_config["simulation_cycles"],
            "NumberOfInitializationCycles"  : simulation_config["initialization_cycles"] if not restart else 0,
            "FrameworkName"                 : material.uuid,
            "HeliumVoidFraction"            : void_fraction,
            "ExternalTemperature"           : simulation_config["temperature"],
            "ExternalPressure"              : simulation_config["pressure"],
            "MoleculeName"                  : simulation_config["adsorbate"],
            "UnitCell"                      : " ".join(map(str, unit_cells))}

    # Load template and replace values
    input_data = load_and_subs_template("input_file_templates/gas_loading.input", values)

    # Write simulation input-file
    with open(filename, "w") as raspa_input_file:
        raspa_input_file.write(input_data)

def write_output_files(material, simulation_config, output_dir, restart=False, filename=None):
    # Write simulation input-files
    # RASPA input-file
    if filename is None:
        filename = os.path.join(output_dir, "{}_loading.input".format(simulation_config['adsorbate']))
    write_raspa_file(filename, material, simulation_config, restart)
    # Pseudomaterial mol-file
    write_mol_file(material, output_dir)
    # Lennard-Jones parameters, force_field_mixing_rules.def
    write_mixing_rules(material.structure, output_dir)
    # Pseudoatom definitions, pseudo_atoms.def (placeholder values)
    write_pseudo_atoms(material.structure, output_dir)
    # Overwritten interactions, force_field.def (none overwritten by default)
    write_force_field(output_dir)

def parse_output(output_file, material, simulation_config):
    """Parse output file for gas adsorption data.

    Args:
        output_file (str): path to simulation output file.

    Returns:
        results (dict): absolute and excess molar, gravimetric, and volumetric
            gas loadings, as well as energy of average, van der Waals, and
            Coulombic host-host, host-adsorbate, and adsorbate-adsorbate
            interactions.

    """
    gas_loading = GasLoading()
    gas_loading.adsorbate        = simulation_config["adsorbate"]
    gas_loading.pressure         = simulation_config["pressure"]
    gas_loading.temperature      = simulation_config["temperature"]
    atom_blocks = []

    with open(output_file) as origin:
        lines = origin.read().split('\n')
        for i, line in enumerate(lines):
            if "absolute [mol/kg" in line:
                gas_loading.absolute_molar_loading = float(line.split()[5])
            elif "absolute [cm^3 (STP)/g" in line:
                gas_loading.absolute_gravimetric_loading = float(line.split()[6])
            elif "absolute [cm^3 (STP)/c" in line:
                gas_loading.absolute_volumetric_loading = float(line.split()[6])
                gas_loading.absolute_volumetric_loading_error = float(line.split()[8])
            elif "excess [mol/kg" in line:
                gas_loading.excess_molar_loading = float(line.split()[5])
            elif "excess [cm^3 (STP)/g" in line:
                gas_loading.excess_gravimetric_loading = float(line.split()[6])
            elif "excess [cm^3 (STP)/c" in line:
                gas_loading.excess_volumetric_loading  = float(line.split()[6])
            elif "Average Host-Host energy:" in line:
                line = lines[i + 8]
                gas_loading.host_host_avg = float(line.split()[1])
                gas_loading.host_host_vdw = float(line.split()[5])
                gas_loading.host_host_cou = float(line.split()[7])
            elif "Average Adsorbate-Adsorbate energy:" in line:
                line = lines[i + 8]
                gas_loading.adsorbate_adsorbate_avg = float(line.split()[1])
                gas_loading.adsorbate_adsorbate_vdw = float(line.split()[5])
                gas_loading.adsorbate_adsorbate_cou = float(line.split()[7])
            elif "Average Host-Adsorbate energy:" in line:
                line = lines[i + 8]
                gas_loading.host_adsorbate_avg  = float(line.split()[1])
                gas_loading.host_adsorbate_vdw  = float(line.split()[5])
                gas_loading.host_adsorbate_cou  = float(line.split()[7])
            elif "Number of molecules:" in line:
                atom_blocks = [float(lines[offset + i + 5].split()[2]) for offset in range(5)]
            elif "Conversion factor molecules/unit cell -> cm^3 STP/cm^3:" in line:
                atoms_uc_to_vv = float(line.split()[7])

        print("{} LOADING : {} v/v (STP)".format(simulation_config["adsorbate"],
                                            gas_loading.absolute_volumetric_loading))
        if material.parent:
            print("(parent LOADING : {} v/v (STP))".format(material.parent.gas_loading[0].absolute_volumetric_loading))


    return gas_loading, atom_blocks, atoms_uc_to_vv

def pressure_string(p):
    if p >= 10 ** 6:
        return "{:.1e}".format(p)
    else:
        return str(p)

def run(material, simulation_config, config):
    """Runs gas loading simulation.

    Args:
        material_id (Material): material record.

    Returns:
        results (dict): gas loading simulation results.

    """
    adsorbate = simulation_config["adsorbate"]
    output_dir = "output_{}_{}".format(material.uuid, uuid4())
    os.makedirs(output_dir, exist_ok=True)
    raspa_config = "./{}_loading.input".format(adsorbate)
    raspa_restart_config = "./{}_loading_restart.input".format(adsorbate)

    # RASPA input-files
    write_output_files(material, simulation_config, output_dir, restart=False, filename=os.path.join(output_dir, raspa_config))
    write_output_files(material, simulation_config, output_dir, restart=True, filename=os.path.join(output_dir, raspa_restart_config))

    # Run simulations
    print("--")
    print("Date             : {}".format(datetime.now().date().isoformat()))
    print("Time             : {}".format(datetime.now().time().isoformat()))
    print("Simulation type  : {}".format(simulation_config["type"]))
    print("Adsorbate        : {}".format(adsorbate))
    print("Pressure         : {}".format(simulation_config["pressure"]))
    print("Temperature      : {}".format(simulation_config["temperature"]))

    unit_cells = material.structure.minimum_unit_cells(simulation_config['cutoff'])
    total_unit_cells = unit_cells[0] * unit_cells[1] * unit_cells[2]
    all_atom_blocks = []

    subprocess.run(["simulate", "-i", raspa_config], check=True, cwd=output_dir)
    for i in range(simulation_config['max_restarts'] + 1):

        data_files = glob(os.path.join(output_dir, "Output", "System_0", "*.data"))
        if len(data_files) != 1:
            raise Exception("ERROR: There should only be one data file in the output directory for %s. Check code!" % output_dir)
        output_file = data_files[0]

        # Parse output
        gas_loading, atom_blocks, atoms_uc_to_vv = parse_output(output_file, material, simulation_config)
        atom_blocks = [a * atoms_uc_to_vv / total_unit_cells for a in atom_blocks]
        print("new blocks for averaging [v/v]: ", atom_blocks)
        print("atoms_uc_to_vv = %f" % atoms_uc_to_vv)
        print("reported V/V: %f" % gas_loading.absolute_volumetric_loading)
        print("reported err: %f" % gas_loading.absolute_volumetric_loading_error)


        all_atom_blocks += atom_blocks
        print("all blocks: ", all_atom_blocks)

        # assign two initialization blocks to every restart run
        run_blocks = all_atom_blocks[math.floor(i/2)*5:]
        print("run blocks: ", run_blocks)
        print("run blocks len: %f" % (len(run_blocks) / 5))

        blocks_for_averaging = np.mean(np.array(run_blocks).reshape(-1, int(len(run_blocks) / 5)), axis=1)
        print("incorporated blocks for averaging [v/v]: ", blocks_for_averaging)
        atoms_std = np.std(blocks_for_averaging)
        print("2*std of all blocks avg %d: %f" % (i, atoms_std*2))
        print("2*std of all blocks: %f" % (2 * np.std(run_blocks)))

        error_vv = 2*atoms_std * atoms_uc_to_vv / total_unit_cells
        gas_loading.absolute_volumetric_loading = np.mean(blocks_for_averaging)
        gas_loading.absolute_volumetric_loading_error = error_vv
        print("calculated V/V: %f" % gas_loading.absolute_volumetric_loading)
        print("calculated error: %f" % error_vv)

        print("Copying restart to RestartInitial...")
        # remove old RestartInitial directory and copy the current one to there
        shutil.rmtree(os.path.join(output_dir, "RestartInitial"), ignore_errors=True)
        shutil.copytree(os.path.join(output_dir, "Restart"), os.path.join(output_dir, "RestartInitial"))

        print("Moving backup RASPA outputs to restart index")
        shutil.move(os.path.join(output_dir, "Output"), os.path.join(output_dir, "Output-%d" % i))
        shutil.move(os.path.join(output_dir, "Restart"), os.path.join(output_dir, "Restart-%d" % i))
        shutil.move(os.path.join(output_dir, "Movies"), os.path.join(output_dir, "Movies-%d" % i))
        shutil.move(os.path.join(output_dir, "VTK"), os.path.join(output_dir, "VTK-%d" % i))

        gas_loading.cycles = simulation_config['simulation_cycles'] * (i + 1)
        if (gas_loading.absolute_volumetric_loading_error < simulation_config['restart_err_threshold']):
            print("Exiting because v/v err < restart_err_threshold: %4.2f < %4.2f" %
                (gas_loading.absolute_volumetric_loading_error, simulation_config['restart_err_threshold']))
            break
        elif i == simulation_config['max_restarts']:
            print("Exiting because we've already restarted maximum number of times.")
            print("v/v err >= restart_err_threshold: %4.2f >= %4.2f" %
                (gas_loading.absolute_volumetric_loading_error, simulation_config['restart_err_threshold']))
            print("--")

            break
        else:
            print("\n--")
            print("restart # %d" % i)
            subprocess.run(["simulate", "-i", raspa_restart_config], check=True, cwd=output_dir)




    material.gas_loading.append(gas_loading)

    if not config['keep_configs']:
        shutil.rmtree(output_dir, ignore_errors=True)
    sys.stdout.flush()
