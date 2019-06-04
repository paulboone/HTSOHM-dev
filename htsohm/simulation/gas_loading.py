from datetime import datetime
from glob import glob
import os
from pathlib import Path
import subprocess
import shutil
from string import Template
import sys
from uuid import uuid4

from htsohm.material_files import write_mol_file, write_mixing_rules
from htsohm.material_files import write_pseudo_atoms, write_force_field
from htsohm.simulation.files import load_and_subs_template
from htsohm.db import GasLoading

def write_raspa_file(filename, material, simulation_config):
    """Writes RASPA input file for simulating gas adsorption.

    Args:
        filename (str): path to input file.
        material (Material): material record.

    Writes RASPA input-file.

    """
    # Set default void fraction value if non was calculated
    if material.void_fraction[0].void_fraction == None:
        void_fraction = 0.
    else:
        void_fraction = material.void_fraction[0].void_fraction

    # Load simulation parameters from config
    unit_cells = material.structure.minimum_unit_cells(simulation_config['cutoff'])
    values = {
            "Cutoff"                        : simulation_config['cutoff'],
            "NumberOfCycles"                : simulation_config["simulation_cycles"],
            "NumberOfInitializationCycles"  : simulation_config["initialization_cycles"],
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

def write_output_files(material, simulation_config, output_dir):
    # Write simulation input-files
    # RASPA input-file
    filename = os.path.join(output_dir, "{}_loading.input".format(simulation_config['adsorbate']))
    write_raspa_file(filename, material, simulation_config)
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

    with open(output_file) as origin:
        line_counter = 1
        for line in origin:
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
                host_host_line = line_counter + 8
            elif "Average Adsorbate-Adsorbate energy:" in line:
                adsorbate_adsorbate_line = line_counter + 8
            elif "Average Host-Adsorbate energy:" in line:
                host_adsorbate_line = line_counter + 8
            line_counter += 1

        print("\n{} LOADING : {} v/v (STP)".format(simulation_config["adsorbate"],
                                            gas_loading.absolute_volumetric_loading))
        if material.parent:
            print("(parent LOADING : {} v/v (STP))".format(material.parent.gas_loading[0].absolute_volumetric_loading))
        print("\n")
    with open(output_file) as origin:
        line_counter = 1
        for line in origin:
            if line_counter == host_host_line:
                gas_loading.host_host_avg = float(line.split()[1])
                gas_loading.host_host_vdw = float(line.split()[5])
                gas_loading.host_host_cou = float(line.split()[7])
            elif line_counter == adsorbate_adsorbate_line:
                gas_loading.adsorbate_adsorbate_avg = float(line.split()[1])
                gas_loading.adsorbate_adsorbate_vdw = float(line.split()[5])
                gas_loading.adsorbate_adsorbate_cou = float(line.split()[7])
            elif line_counter == host_adsorbate_line:
                gas_loading.host_adsorbate_avg  = float(line.split()[1])
                gas_loading.host_adsorbate_vdw  = float(line.split()[5])
                gas_loading.host_adsorbate_cou  = float(line.split()[7])
            line_counter += 1

    material.gas_loading.append(gas_loading)

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

    # RASPA input-file
    write_output_files(material, simulation_config, output_dir)

    # Run simulations
    print("Date             : {}".format(datetime.now().date().isoformat()))
    print("Time             : {}".format(datetime.now().time().isoformat()))
    print("Simulation type  : {}".format(simulation_config["type"]))
    print("Adsorbate        : {}".format(adsorbate))
    print("Pressure         : {}".format(simulation_config["pressure"]))
    print("Temperature      : {}".format(simulation_config["temperature"]))

    subprocess.run(["simulate", "-i", "./{}_loading.input".format(adsorbate)], check=True, cwd=output_dir)

    data_files = glob(os.path.join(output_dir, "Output", "System_0", "*.data"))
    if len(data_files) != 1:
        raise Exception("ERROR: There should only be one data file in the output directory for %s. Check code!" % output_dir)
    output_file = data_files[0]

    # Parse output
    parse_output(output_file, material, simulation_config)
    if not config['keep_configs']:
        shutil.rmtree(output_dir, ignore_errors=True)
    sys.stdout.flush()
