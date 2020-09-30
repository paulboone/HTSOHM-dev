from glob import glob
import os
import re
import subprocess
import shutil
from string import Template
import sys

import numpy as np

from htsohm.simulation.raspa import write_mol_file, write_mixing_rules, write_pseudo_atoms, write_force_field
from htsohm.simulation.templates import load_and_subs_template
from htsohm.db import HenrysCoefficient
from htsohm.slog import slog

def write_raspa_file(filename, material, simulation_config, restart):
    # Load simulation parameters from config
    unit_cells = material.structure.minimum_unit_cells(simulation_config['cutoff'])
    values = {
            "Restart"                       : 'yes' if restart else 'no',
            "Cutoff"                        : simulation_config['cutoff'],
            "NumberOfCycles"                : simulation_config["simulation_cycles"],
            "NumberOfInitializationCycles"  : simulation_config["initialization_cycles"] if not restart else 0,
            "FrameworkName"                 : material.id,
            "ExternalTemperature"           : simulation_config["temperature"],
            "UnitCell"                      : " ".join(map(str, unit_cells)),
            "RosenbluthWeight"              : simulation_config["rosenbluth_weight"]}
            # "MoleculeName"                  : simulation_config["adsorbate"],

    # Load template and replace values
    input_data = load_and_subs_template("input_file_templates/henrys_coefficients.input", values)

    # Write simulation input-file
    with open(filename, "w") as raspa_input_file:
        raspa_input_file.write(input_data)

def write_input_files(material, simulation_config, output_dir, restart=False, filename=None):
    # Write simulation input-files
    # RASPA input-file
    if filename is None:
        filename = os.path.join(output_dir, "henrys.input")
    write_raspa_file(filename, material, simulation_config, restart)
    # Pseudomaterial mol-file
    write_mol_file(material, output_dir)
    # Lennard-Jones parameters, force_field_mixing_rules.def
    write_mixing_rules(material.structure, output_dir)
    # Pseudoatom definitions, pseudo_atoms.def (placeholder values)
    write_pseudo_atoms(material.structure, output_dir)
    # Overwritten interactions, force_field.def (none overwritten by default)
    write_force_field(output_dir)

def parse_output(output_file, simulation_config):
    fl = r"[-+]?\d+(?:\.\d*)?(?:[eE][-+]?\d+)?"

    # captured groups are [gas, Henry's coefficient, Henry's coefficient error]
    henrys_re = re.compile(r"\s+\[(\w+)\] Average Henry coefficient:  ({fl}) \+\/\- ({fl}) \[mol/kg/Pa\]".format(fl=fl))

    gas_henrys_error = []
    with open(output_file) as f:
        for line in f:
            m=re.match(henrys_re, line)
            if m:
                gas_henrys_error += [m.groups()]
    return gas_henrys_error

def run(material, simulation_config, config):
    output_dir = "output_{}_{}".format(material.id, simulation_config['name'])
    os.makedirs(output_dir, exist_ok=True)
    slog("Output directory : {}".format(output_dir))
    raspa_config = "./henrys_coefficients.input"

    # RASPA input-files
    write_input_files(material, simulation_config, output_dir, restart=False, filename=os.path.join(output_dir, raspa_config))

    # Run simulations
    process = subprocess.run(["simulate", "-i", raspa_config], check=True, cwd=output_dir, capture_output=True, text=True)
    slog(process.stdout)

    data_files = glob(os.path.join(output_dir, "Output", "System_0", "*.data"))
    if len(data_files) != 1:
        raise Exception("ERROR: There should only be one data file in the output directory for %s. Check code!" % output_dir)

    gas_henrys_error = parse_output(data_files[0], simulation_config)


    hc = HenrysCoefficient()
    for (gas, henrys, henrys_error) in gas_henrys_error:
        if gas == "CO2":
            hc.co2_henrys = henrys
            hc.co2_henrys_error = henrys_error
        elif gas == "N2":
            hc.n2_henrys = henrys
            hc.n2_henrys_error = henrys_error
        elif gas == "water":
            hc.h2o_henrys = henrys
            hc.h2o_henrys_error = henrys_error
    material.henrys_coefficient.append(hc)

    if not config['keep_configs']:
        shutil.rmtree(output_dir, ignore_errors=True)
    sys.stdout.flush()
