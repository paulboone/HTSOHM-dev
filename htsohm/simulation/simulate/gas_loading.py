from datetime import datetime
from glob import glob
import math
import os
from pathlib import Path
import subprocess
import shutil
from string import Template
import sys

import numpy as np

from htsohm.simulation.raspa import write_mol_file, write_mixing_rules
from htsohm.simulation.raspa import write_pseudo_atoms, write_force_field
from htsohm.simulation.templates import load_and_subs_template
from htsohm.slog import slog

def write_raspa_file(filename, material, simulation_config, restart):
    # Set default void fraction value if non was calculated
    if len(material.void_fraction) == 0 or material.void_fraction[0].void_fraction is None:
        void_fraction = 0.
    else:
        void_fraction = material.void_fraction[0].void_fraction

    # Load simulation parameters from config
    unit_cells = material.minimum_unit_cells(simulation_config['cutoff'])
    values = {
            "Restart"                       : 'yes' if restart else 'no',
            "Cutoff"                        : simulation_config['cutoff'],
            "NumberOfCycles"                : simulation_config["simulation_cycles"],
            "NumberOfInitializationCycles"  : simulation_config["initialization_cycles"] if not restart else 0,
            "FrameworkName"                 : material.id_or_uuid,
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

def write_input_files(material, simulation_config, output_dir, restart=False, filename=None):
    # Write simulation input-files
    # RASPA input-file
    if filename is None:
        filename = os.path.join(output_dir, "{}_loading.input".format(simulation_config['adsorbate']))
    write_raspa_file(filename, material, simulation_config, restart)
    # Pseudomaterial mol-file
    write_mol_file(material, output_dir)
    # Lennard-Jones parameters, force_field_mixing_rules.def
    write_mixing_rules(material, output_dir)
    # Pseudoatom definitions, pseudo_atoms.def (placeholder values)
    write_pseudo_atoms(material, output_dir)
    # Overwritten interactions, force_field.def (none overwritten by default)
    write_force_field(output_dir)

def parse_output(output_file, material, simulation_config):
    gas_loading = GasLoading()
    gas_loading.adsorbate        = simulation_config["adsorbate"]
    gas_loading.pressure         = simulation_config["pressure"]
    gas_loading.temperature      = simulation_config["temperature"]
    atom_blocks = []

    with open(output_file) as origin:
        lines = origin.read().split('\n')
        for i, line in enumerate(lines):

            if "absolute [cm^3 (STP)/c" in line:
                gas_loading.absolute_volumetric_loading = float(line.split()[6])
                gas_loading.absolute_volumetric_loading_error = float(line.split()[8])
            elif "Number of molecules:" in line:
                atom_blocks = [float(lines[offset + i + 5].split()[2]) for offset in range(5)]
            elif "Conversion factor molecules/unit cell -> cm^3 STP/cm^3:" in line:
                atoms_uc_to_vv = float(line.split()[7])

        slog("{} LOADING : {} v/v (STP)".format(simulation_config["adsorbate"],
                                            gas_loading.absolute_volumetric_loading))
        if material.parent:
            slog("(parent LOADING : {} v/v (STP))".format(material.parent.gas_loading[0].absolute_volumetric_loading))


    return gas_loading, atom_blocks, atoms_uc_to_vv

def pressure_string(p):
    if p >= 10 ** 6:
        return "{:.1e}".format(p)
    else:
        return str(p)

def run(material, simulation_config, config):
    raise "need updating to new database format"
    adsorbate = simulation_config["adsorbate"]
    output_dir = "output_{}_{}".format(material.id_or_uuid, simulation_config['name'])
    os.makedirs(output_dir, exist_ok=True)
    raspa_config = "./{}_loading.input".format(adsorbate)
    raspa_restart_config = "./{}_loading_restart.input".format(adsorbate)

    # RASPA input-files
    write_input_files(material, simulation_config, output_dir, restart=False, filename=os.path.join(output_dir, raspa_config))
    write_input_files(material, simulation_config, output_dir, restart=True, filename=os.path.join(output_dir, raspa_restart_config))

    # Run simulations
    slog("Adsorbate        : {}".format(adsorbate))
    slog("Pressure         : {}".format(simulation_config["pressure"]))
    slog("Temperature      : {}".format(simulation_config["temperature"]))

    unit_cells = material.minimum_unit_cells(simulation_config['cutoff'])
    total_unit_cells = unit_cells[0] * unit_cells[1] * unit_cells[2]
    all_atom_blocks = []

    process = subprocess.run(["simulate", "-i", raspa_config], cwd=output_dir, capture_output=True, text=True)
    if process.returncode != 0:
        slog(process.stdout)
        slog(process.stderr)
        process.check_returncode()

    for i in range(simulation_config['max_restarts'] + 1):

        data_files = glob(os.path.join(output_dir, "Output", "System_0", "*.data"))
        if len(data_files) != 1:
            raise Exception("ERROR: There should only be one data file in the output directory for %s. Check code!" % output_dir)
        output_file = data_files[0]

        # Parse output
        gas_loading, atom_blocks, atoms_uc_to_vv = parse_output(output_file, material, simulation_config)
        atom_blocks = [a * atoms_uc_to_vv / total_unit_cells for a in atom_blocks]
        slog("new blocks for averaging [v/v]: ", atom_blocks)
        slog("atoms_uc_to_vv = %f" % atoms_uc_to_vv)
        slog("reported V/V: %f" % gas_loading.absolute_volumetric_loading)
        slog("reported err: %f" % gas_loading.absolute_volumetric_loading_error)


        all_atom_blocks += atom_blocks
        slog("all blocks: ", all_atom_blocks)

        # assign two initialization blocks to every restart run
        run_blocks = all_atom_blocks[math.floor(i/2)*5:]
        slog("run blocks: ", run_blocks)
        slog("run blocks len: %f" % (len(run_blocks) / 5))

        blocks_for_averaging = np.mean(np.array(run_blocks).reshape(-1, int(len(run_blocks) / 5)), axis=1)
        slog("incorporated blocks for averaging [v/v]: ", blocks_for_averaging)
        atoms_std = np.std(blocks_for_averaging)
        slog("2*std of all blocks avg %d: %f" % (i, atoms_std*2))
        slog("2*std of all blocks: %f" % (2 * np.std(run_blocks)))

        error_vv = 2*atoms_std * atoms_uc_to_vv / total_unit_cells
        gas_loading.absolute_volumetric_loading = np.mean(blocks_for_averaging)
        gas_loading.absolute_volumetric_loading_error = error_vv
        slog("calculated V/V: %f" % gas_loading.absolute_volumetric_loading)
        slog("calculated error: %f" % error_vv)

        slog("Copying restart to RestartInitial...")
        # remove old RestartInitial directory and copy the current one to there
        shutil.rmtree(os.path.join(output_dir, "RestartInitial"), ignore_errors=True)
        shutil.copytree(os.path.join(output_dir, "Restart"), os.path.join(output_dir, "RestartInitial"))

        slog("Moving backup RASPA outputs to restart index")
        shutil.move(os.path.join(output_dir, "Output"), os.path.join(output_dir, "Output-%d" % i))
        shutil.move(os.path.join(output_dir, "Restart"), os.path.join(output_dir, "Restart-%d" % i))
        shutil.move(os.path.join(output_dir, "Movies"), os.path.join(output_dir, "Movies-%d" % i))
        shutil.move(os.path.join(output_dir, "VTK"), os.path.join(output_dir, "VTK-%d" % i))

        gas_loading.cycles = simulation_config['simulation_cycles'] * (i + 1)
        if (gas_loading.absolute_volumetric_loading_error < simulation_config['restart_err_threshold']):
            slog("Exiting because v/v err < restart_err_threshold: %4.2f < %4.2f" %
                (gas_loading.absolute_volumetric_loading_error, simulation_config['restart_err_threshold']))
            break
        elif i == simulation_config['max_restarts']:
            slog("Exiting because we've already restarted maximum number of times.")
            slog("v/v err >= restart_err_threshold: %4.2f >= %4.2f" %
                (gas_loading.absolute_volumetric_loading_error, simulation_config['restart_err_threshold']))
            break
        else:
            slog("\n--")
            slog("restart # %d" % i)
            process = subprocess.run(["simulate", "-i", raspa_restart_config], cwd=output_dir, capture_output=True, text=True)
            if process.returncode != 0:
                slog(process.stdout)
                slog(process.stderr)
                process.check_returncode()




    material.gas_loading.append(gas_loading)

    if not config['keep_configs']:
        shutil.rmtree(output_dir, ignore_errors=True)
    sys.stdout.flush()
