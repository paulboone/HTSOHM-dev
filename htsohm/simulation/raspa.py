import os

def write_mol_file(material, simulation_path):
    """Writes .mol file for structural information."""

    file_name = os.path.join(simulation_path, "{}.mol".format(material.id_or_uuid))
    with open(file_name, "w") as mol_file:
        mol_file.write(
                " Molecule_name: {}\n".format(material.id_or_uuid) +
                "\n" +
                "  Coord_Info: Listed Cartesian None\n" +
                "        {}\n".format(len(material.atom_sites)))
        for i in range(len(material.atom_sites)):
            a = material.atom_sites[i]
            mol_file.write(
                    "{:6} {:10.4f} {:10.4f} {:10.4f}  {:5} {:10.8f}  0  0\n".format(
                        i + 1, round(a.x * material.a, 4),
                        round(a.y * material.b, 4),
                        round(a.z * material.c, 4),
                        str(a.atom_types.atom_type_index()), round(a.q, 8)))
        mol_file.write(
                "\n" +
                "\n" +
                "\n" +
                "  Fundcell_Info: Listed\n" +
                "        {:10.4f}       {:10.4f}       {:10.4f}\n".format(
                    round(material.a, 4), round(material.b, 4), round(material.c, 4)) +
                "           90.0000          90.0000          90.0000\n" +
                "           0.00000          0.00000          0.00000\n" +
                "        {:10.4f}       {:10.4f}       {:10.4f}\n".format(
                    round(material.a, 4), round(material.b, 4), round(material.c, 4)) +
                "\n")

def write_mixing_rules(material, simulation_path):
    """Writes .def file for forcefield information."""
    adsorbate_LJ_atoms = [
            ['N_n2',    36.0,       3.31],
            ['C_co2',   27.0,       2.80],
            ['O_co2',   79.0,       3.05],
            ['CH4_sp3', 158.5,      3.72],
            ['He',      10.9,       2.64],
            ['H_com',   36.7,       2.958],
            ['Kr',      167.06,     3.924],
            ['Xe',      110.704,    3.690],
            ['Ow',      78.0,     3.154]
    ]

    adsorbate_none_atoms = ['N_com', 'H_h2', 'Lw', 'Hw']

    file_name = os.path.join(simulation_path, 'force_field_mixing_rules.def')
    with open(file_name, "w") as f:
        f.write("""# general rule for shifted vs truncated
shifted
# general rule tailcorrections
no
# number of defined interactions
{0}
# type interaction, parameters\n""".format(len(material.atom_sites) + len(adsorbate_LJ_atoms) + len(adsorbate_none_atoms)))

        # write one atom type per atom site, so we can define per-site charges on the types.
        type_template = "{0:12} lennard-jones {1:.4f} {2:.4f}\n"
        for i, a in enumerate(material.atom_sites):
            f.write(type_template.format(i, a.atom_types.epsilon, a.atom_types.sigma))
        for at in adsorbate_LJ_atoms:
            f.write(type_template.format(at[0], at[1], at[2]))

        for at in adsorbate_none_atoms:
            f.write(
                "{0:12} none\n".format(at)
            )
        f.write(
            "# general mixing rule for Lennard-Jones\n" +
            "Lorentz-Berthelot")

def write_pseudo_atoms(material, simulation_path):
    """Writes .def file for chemical information.

    Args:
        simulation_path (str): path to pseudo atoms definitions file

    Returns:
        NOTE: ALL CHARGES ARE 0. IN THIS VERSION.

    """
    temporary_charge = 0.

    file_name = os.path.join(simulation_path, 'pseudo_atoms.def')
    with open(file_name, "w") as pseudo_atoms_file:
        pseudo_atoms_file.write(
            "#number of pseudo atoms\n" +
            "%s\n" % (len(material.atom_sites) + 13) +
            "#type  print   as  chem    oxidation   mass    charge  polarization    B-factor    radii   " +
                 "connectivity  anisotropic anisotrop-type  tinker-type\n")

        atom_type_string = "{:7}  yes  C   C   0   12.0       {:f}  0.0  1.0  1.0    0  0  absolute  0\n"
        for i, a in enumerate(material.atom_sites):
            pseudo_atoms_file.write(atom_type_string.format(i, a.q))

        pseudo_atoms_file.write(
            "N_n2     yes  N   N   0   14.00674   -0.4048   0.0  1.0  0.7    0  0  relative  0\n" +
            "N_com    no   N   -   0    0.0        0.8096   0.0  1.0  0.7    0  0  relative  0\n" +
            "C_co2    yes  C   C   0   12.0        0.70     0.0  1.0  0.720  0  0  relative  0\n" +
            "O_co2    yes  O   O   0   15.9994    -0.35     0.0  1.0  0.68   0  0  relative  0\n" +
            "CH4_sp3  yes  C   C   0   16.04246    0.0      0.0  1.0  1.00   0  0  relative  0\n" +
            "He       yes  He  He  0    4.002602   0.0      0.0  1.0  1.0    0  0  relative  0\n" +
            "H_h2     yes  H   H   0    1.00794    0.468    0.0  1.0  0.7    0  0  relative  0\n" +
            "H_com    no   H   H   0    0.0        0.936    0.0  1.0  0.7    0  0  relative  0\n" +
            "Xe       yes  Xe  Xe  0  131.293      0.0      0.0  1.0  2.459  0  0  relative  0\n" +
            "Kr       yes  Kr  Kr  0   83.798      0.0      0.0  1.0  2.27   0  0  relative  0\n" +
            "Ow       yes  O   O   0   15.9996     0.0      0.0  1.0  0.5    2  0  absolute  0\n" +
            "Hw       yes  H   H   0    1.0008     0.52     0.0  1.0  1.00   1  0  absolute  0\n" +
            "Lw       no   L   -   0    0.0       -1.04     0.0  1.0  1.00   1  0  absolute  0\n"
        )

def write_force_field(simulation_path):
    """Writes .def file to overwrite LJ-type interactions.

    Args:
        file_name (str): path to write .def-file

    NOTE: NO INTERACTIONS ARE OVERWRITTEN BY DEFAULT.

    """
    forcefield_file = """# rules to overwrite
0
# number of defined interactions
0
# mixing rules to overwrite
0"""

    file_name = os.path.join(simulation_path, 'force_field.def')
    with open(file_name, "w") as f:
        f.write(forcefield_file)
