
from glob import glob
import math
import os
import shutil
import sys
import subprocess
import time

from datetime import datetime
from string import Template
from pathlib import Path

from htsohm.simulation.raspa import write_mol_file, write_mixing_rules
from htsohm.simulation.raspa import write_pseudo_atoms, write_force_field
from htsohm.simulation.templates import load_and_subs_template
from htsohm.void_fraction import calculate_void_fraction
from htsohm.slog import slog

def write_raspa_file(filename, material, simulation_config):
    # Load simulation parameters from config
    unit_cells = material.minimum_unit_cells(simulation_config['cutoff'])
    values = {
            "Cutoff"                 : simulation_config['cutoff'],
            "NumberOfCycles"         : simulation_config["simulation_cycles"],
            "FrameworkName"          : material.id_or_uuid,
            "ExternalTemperature"    : simulation_config["temperature"],
            "MoleculeName"           : simulation_config["adsorbate"],
            "UnitCell"               : " ".join(map(str, unit_cells))}

    # Load template and replace values
    input_data = load_and_subs_template("input_file_templates/void_fraction.input", values)

    # Write simulation input-file
    with open(filename, "w") as raspa_input_file:
        raspa_input_file.write(input_data)

def write_input_files(material, simulation_config, output_dir):
    # Write simulation input-files
    # RASPA input-file
    filename = os.path.join(output_dir, "void_fraction.input")
    write_raspa_file(filename, material, simulation_config)
    # Pseudomaterial mol-file
    write_mol_file(material, output_dir)
    # Lennard-Jones parameters, force_field_mixing_rules.def
    write_mixing_rules(material, output_dir)
    # Pseudoatom definitions, pseudo_atoms.def (placeholder values)
    write_pseudo_atoms(material, output_dir)
    # Overwritten interactions, force_field.def (none overwritten by default)
    write_force_field(output_dir)

def parse_output(output_file):
    vf = None
    with open(output_file) as origin:
        for line in origin:
            if not "Average Widom Rosenbluth-weight:" in line:
                continue
            vf = float(line.split()[4])
    return vf


def run(material, simulation_config, config):
    output_dir = "output_{}{}".format(simulation_config['prefix'], material.id_or_uuid)
    slog("Output directory : {}".format(output_dir))
    os.makedirs(output_dir, exist_ok=True)

    if "raspa" in simulation_config["fields"]:
        slog("Temperature      : {}".format(simulation_config["raspa"]["temperature"]))
        slog("Probe            : {}".format(simulation_config["raspa"]["adsorbate"]))

        write_input_files(material, simulation_config["raspa"], output_dir)
        tbegin = time.perf_counter()
        process = subprocess.run(["simulate", "-i", "./void_fraction.input"], cwd=output_dir, capture_output=True, text=True)
        if process.returncode != 0:
            slog(process.stdout)
            slog(process.stderr)
            process.check_returncode()

        data_files = glob(os.path.join(output_dir, "Output", "System_0", "*.data"))
        if len(data_files) != 1:
            raise Exception("ERROR: There should only be one data file in the output directory for %s. Check code!" % output_dir)
        output_file = data_files[0]

        vf_raspa = parse_output(output_file)
        material.__setattr__("%sraspa" % simulation_config["prefix"], vf_raspa)
        slog("RASPA void fraction simulation time: %5.2f seconds" % (time.perf_counter() - tbegin))
        slog("RASPA VOID FRACTION : {}".format(vf_raspa))

    # run geometric void fraction
    if "geo" in simulation_config["fields"]:
        slog("Probe radius [geo]: {}".format(simulation_config["geo"]["probe_radius"]))
        tbegin = time.perf_counter()
        atoms = [(a.x * material.a, a.y * material.b, a.z * material.c, a.atom_types.sigma) for a in material.atom_sites]
        box = (material.a, material.b, material.c)
        vf_geo = calculate_void_fraction(atoms, box, probe_r=simulation_config["geo"]["probe_radius"])
        material.__setattr__("%sgeo" % simulation_config["prefix"], vf_geo)
        slog("GEOMETRIC void fraction: %f" % vf_geo)
        slog("GEOMETRIC void fraction simulation time: %5.2f   seconds" % (time.perf_counter() - tbegin))

    if not config['keep_configs']:
        shutil.rmtree(output_dir, ignore_errors=True)
    sys.stdout.flush()
