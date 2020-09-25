import sys
import os
import subprocess
import shutil
from datetime import datetime
from uuid import uuid4
from string import Template
from pathlib import Path

from htsohm.simulation.raspa import write_mol_file, write_mixing_rules
from htsohm.simulation.raspa import write_pseudo_atoms, write_force_field
from htsohm.simulation.templates import load_and_subs_template
from htsohm.db import SurfaceArea
from htsohm.slog import slog

def write_raspa_file(filename, material, simulation_config):
    """Writes RASPA input file for calculating surface area.

    Args:
        filename (str): path to input file.
        material_id (str): uuid for material.

    Writes RASPA input-file.

    """
    # Load simulation parameters from config
    unit_cells = material.structure.minimum_unit_cells(simulation_config['cutoff'])
    values = {
            "Cutoff"            : simulation_config['cutoff'],
            "NumberOfCycles"    : simulation_config["simulation_cycles"],
            "FrameworkName"     : material.uuid,
            "MoleculeName"      : simulation_config["adsorbate"],
            "UnitCell"          : " ".join(map(str, unit_cells))}

    # Load template and replace values
    input_data = load_and_subs_template("input_file_templates/surface_area.input", values)

    # Write simulation input-file
    with open(filename, "w") as raspa_input_file:
        raspa_input_file.write(input_data)

def write_input_files(material, simulation_config, output_dir):
    # Write simulation input-files
    # RASPA input-file
    filename = os.path.join(output_dir, "SurfaceArea.input")
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
    """Parse output file for void fraction data.

    Args:
        output_file (str): path to simulation output file.

    Returns:
        results (dict): total unit cell, gravimetric, and volumetric surface
            areas.

    """
    surface_area = SurfaceArea()
    surface_area.adsorbate = simulation_config["adsorbate"]

    with open(output_file) as origin:
        count = 0
        for line in origin:
            if "Surface area" in line:
                if count == 0:
                    surface_area.unit_cell_surface_area = float(line.split()[2])
                    count = count + 1
                elif count == 1:
                    surface_area.gravimetric_surface_area = float(line.split()[2])
                    count = count + 1
                elif count == 2:
                    surface_area.volumetric_surface_area = float(line.split()[2])

    slog("\nSURFACE AREA : {} m^2/cm^3\n".format(surface_area.volumetric_surface_area))

    material.surface_area.append(surface_area)

def run(material, simulation_config, config):
    """Runs surface area simulation.

    Args:
        material (Material): material record.

    Returns:
        results (dict): surface area simulation results.

    """
    output_dir = "output_{}_{}".format(material.uuid, uuid4())
    slog("Output directory :\t{}".format(output_dir))
    os.makedirs(output_dir, exist_ok=True)

    # Write simulation input-files
    write_input_files(filename, material, simulation_config)

    # Run simulations
    while True:
        try:
            slog("Probe            : {}".format(simulation_config["adsorbate"]))
            filename = "output_{}_2.2.2_298.000000_0.data".format(material.uuid)
            output_file = os.path.join(output_dir, "Output", "System_0", filename)

            while not Path(output_file).exists():
                process = subprocess.run(["simulate", "-i", "./SurfaceArea.input"], check=True,
                        cwd=output_dir, capture_output=True, text=True)
                slog(process.stdout)

            # Parse output
            parse_output(output_file, material, simulation_config)
            if not config['keep_configs']:
                shutil.rmtree(output_dir, ignore_errors=True)
            sys.stdout.flush()
        except (FileNotFoundError, KeyError) as err:
            slog(err)
            slog(err.args)
            continue
        break
