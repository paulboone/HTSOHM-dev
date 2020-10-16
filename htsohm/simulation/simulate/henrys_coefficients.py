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
from htsohm.slog import slog

def write_raspa_file(filename, material, simulation_config, restart):
    # Load simulation parameters from config
    unit_cells = material.minimum_unit_cells(simulation_config['cutoff'])
    values = {
            "Restart"                       : 'yes' if restart else 'no',
            "Cutoff"                        : simulation_config['cutoff'],
            "NumberOfCycles"                : simulation_config["simulation_cycles"],
            "NumberOfInitializationCycles"  : simulation_config["initialization_cycles"] if not restart else 0,
            "FrameworkName"                 : material.id_or_uuid,
            "ExternalTemperature"           : simulation_config["temperature"],
            "UnitCell"                      : " ".join(map(str, unit_cells)),
            "RosenbluthWeight"              : simulation_config["rosenbluth_weights"]}

    # Load template and replace values
    input_data = load_and_subs_template("input_file_templates/henrys_coefficients.input", values)

    for i, ads in enumerate(simulation_config["adsorbates"]):
        input_data += """
Component %i MoleculeName  %s
MoleculeDefinition        TraPPE
IdealGasRosenbluthWeight  %s
WidomProbability          1.0
CreateNumberOfMolecules   0
""" % (i, ads, simulation_config["rosenbluth_weights"][ads])

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
    write_mixing_rules(material, output_dir)
    # Pseudoatom definitions, pseudo_atoms.def (placeholder values)
    write_pseudo_atoms(material, output_dir)
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
                matches = list(m.groups())
                matches[1] = float(matches[1])
                matches[2] = float(matches[2])
                gas_henrys_error += [matches]
    return gas_henrys_error

def henrys_m_to_v(volume, num_atom_sites, atom_site_mass=12.0):
    """
        henrys_mass in  mol / (kg * Pa)
        volume in cubic angstroms
        atom_site_mass in g
    """
    na = 6.02E+23
    return (num_atom_sites / na) * (atom_site_mass / 1000) / (volume * 1e-30)


def run(material, simulation_config, config):
    output_dir = "output_{}{}".format(simulation_config['prefix'], material.id_or_uuid)
    os.makedirs(output_dir, exist_ok=True)
    slog("Output directory : {}".format(output_dir))
    raspa_config = "./henrys_coefficients.input"

    # RASPA input-files
    write_input_files(material, simulation_config, output_dir, restart=False, filename=os.path.join(output_dir, raspa_config))

    # Run simulations
    process = subprocess.run(["simulate", "-i", raspa_config], cwd=output_dir, capture_output=True, text=True)
    if process.returncode != 0:
        slog(process.stdout)
        slog(process.stderr)
        process.check_returncode()

    data_files = glob(os.path.join(output_dir, "Output", "System_0", "*.data"))
    if len(data_files) != 1:
        raise Exception("ERROR: There should only be one data file in the output directory for %s. Check code!" % output_dir)


    m_to_v = henrys_m_to_v(material.volume, len(material.atom_sites))
    gas_henrys_error = parse_output(data_files[0], simulation_config)

    for gas, henrys, henrys_error in gas_henrys_error:
        if gas in simulation_config['adsorbates']:
            if gas in simulation_config["fields"]:
                material.__setattr__("%s%s" % (simulation_config["prefix"], gas), henrys)
            if ("%sv" % gas) in simulation_config["fields"]:
                material.__setattr__("%s%sv" % (simulation_config["prefix"], gas), henrys * m_to_v)
            if ("%serror" % gas) in simulation_config["fields"]:
                material.__setattr__("%s%serror" % (simulation_config["prefix"], gas), henrys_error)
            if ("%sverror" % gas) in simulation_config["fields"]:
                material.__setattr__("%s%sverror" % (simulation_config["prefix"], gas), henrys_error * m_to_v)

    if not config['keep_configs']:
        shutil.rmtree(output_dir, ignore_errors=True)
    sys.stdout.flush()
