import os
import shutil
import yaml

import numpy as np
from random import random, choice

from htsohm import binning as bng
from htsohm import simulate as sim
from htsohm.runDB_declarative import session, RunData

def first_s(run_id, strength_0):       # Creates `strength` array
    bins = bng.check_number_of_bins(run_id)
    s_array = np.zeros([bins, bins, bins])

    for i in range(bins):
        for j in range(bins):
            for k in range(bins):
                s_array[i, j, k] = strength_0

    wd = os.environ['HTSOHM_DIR']
    np.save(wd + '/config/' + run_id, s_array)

def calculate_s(run_id, generation):

    wd = os.environ['HTSOHM_DIR']

    with open(wd + '/config/' + run_id + '.yaml') as yaml_file:
        config = yaml.load(yaml_file)
    children_per_generation = config["children-per-generation"]
    first = generation * children_per_generation
    last = (generation + 1) * children_per_generation
    child_ids = np.arange(first, last)
    seed_ids = np.arange(0, children_per_generation)

    strength_0 = np.load(wd + '/config/' + run_id + '.npy')       # Load strength-parameter array

    parent_list = []
    for i in child_ids:
        child = session.query(RunData).filter(RunData.run_id == run_id, 
            RunData.material_id == str(i))
        for item in child:
            parent_id = item.parent_id
        parent = session.query(RunData).get(parent_id)
        parent_ml_bin = parent.methane_loading_bin
        parent_sa_bin = parent.surface_area_bin
        parent_vf_bin = parent.void_fraction_bin
        parent_bin = [parent_ml_bin, parent_sa_bin, parent_vf_bin]
        parent_data = [parent_id, parent_bin]
        parent_list.append(parent_data)

    parent_generation_ids = np.arange(first - children_per_generation, 
        last - children_per_generation)
    grandparent_list = []
    all_parents_list = []                        # this isn't really `all` parents, just those in the same bin as a selected parent
    for i in parent_generation_ids:
        parent = session.query(RunData).filter(RunData.run_id == run_id,
            RunData.material_id == str(i))
        for item in parent:
            grandparent_id = item.parent_id
        grandparent = session.query(RunData).get(grandparent_id)
        grandparent_ml_bin = grandparent.methane_loading_bin
        grandparent_sa_bin = grandparent.surface_area_bin
        grandparent_vf_bin = grandparent.void_fraction_bin
        grandparent_bin = [grandparent_ml_bin, grandparent_sa_bin,
            grandparent_vf_bin]
        grandparent_data = [grandparent_id, grandparent_bin]
        grandparent_list.append(grandparent_data)
    
    grandparent_list_redundant = []
    for i in parent_list:
        for j in grandparent_list:
            if i[1] == j[1]:
                data = [i, j]
                grandparent_list_redundant.append(data)

    grandparent_list_clean = []
    for i in grandparent_list_redundant:
        if i not in grandparent_list_clean:
            grandparent_list_clean.append(i)

    counts = bng.count_all(run_id)

    bin_counts = []
    for i in grandparent_list_clean:
        parent_id = i[1][0]                         # here `parent` refers to the parent of a selected parent
        parent_bin = i[1][1]
        child_generation = generation - 1
        parent_bin_count = session.query(RunData).filter(
            RunData.run_id == run_id, RunData.generation == child_generation,
            RunData.methane_loading_bin == parent_bin[0],
            RunData.surface_area_bin == parent_bin[1],
            RunData.void_fraction_bin == parent_bin[2]).count()
        children = session.query(RunData).filter(
            RunData.run_id == run_id, RunData.generation == child_generation,
            RunData.parent_id == parent_id).all()
        children_bins = []
        children_bin_counts = []
        for material in children:
            child_bin = [
                material.methane_loading_bin,
                material.surface_area_bin,
                material.void_fraction_bin
                        ]
            if child_bin != parent_bin:
                children_bins.append(child_bin)
                child_bin_count = session.query(RunData).filter(
                     RunData.run_id == run_id,
                     RunData.generation == child_generation,
                     RunData.methane_loading_bin == child_bin[0],
                     RunData.surface_area_bin == child_bin[1],
                     RunData.void_fraction_bin == child_bin[2]).count()
                children_bin_counts.append(child_bin_count)

        data = [[parent_bin, parent_bin_count], [children_bins, 
        children_bin_counts]]
        bin_counts.append(data)

    for i in bin_counts:
        parent_bin = i[0][0]
        parent_count = i[0][1]
        child_bins = i[1][0]
        child_counts = i[1][1]
        a = parent_bin[0]
        b = parent_bin[1]
        c = parent_bin[2]
        s_0 = strength_0[a, b, c]
        if parent_bin not in child_bins:
            strength_0[a, b, c] = 0.5 * s_0
        if parent_bin in child_bins:
            if parent_count < 1.1 * min(child_counts):
                strength_0[a, b, c] = 0.5 * s_0
            pos = child_bins.index(parent_bin)
            val = 0
            for j in range(len(child_bins)):
                if j != pos:
                    if child_counts[j] > val:
                        val = child_counts[j]
            if parent_count >= 3 * val:
                strength_0[a, b, c] = 1.5 * s_0

    strength_array_file = "%s/config/%s.npy" % (wd, run_id)
    os.remove(strength_array_file)
    np.save(strength_array_file, strength_0)

# function for finding "closest" distance over periodic boundaries
def closest_distance(x_o, x_r):
    a = 1 - x_r + x_o
    b = abs(x_r - x_o)
    c = 1 - x_o + x_r
    dx = min(a, b, c)

    return dx # is closest in R1 the closest in R3???

# given intial and "random" x-fraction, returns new x fraction
def delta_x(x_o, x_r, strength):  # removed random()
    dx = closest_distance(x_o, x_r)
    if (x_o > x_r
            and (x_o - x_r) > 0.5):
        xfrac = round((x_o + strength * dx) % 1., 4)
    if (x_o < x_r
            and (x_r - x_o) > 0.5):
        xfrac = round((x_o - strength * dx) % 1., 4)
    if (x_o >= x_r
            and (x_o - x_r) < 0.5):
        xfrac = round(x_o - strength * dx, 4)
    if (x_o < x_r
            and (x_r - x_o) < 0.5):
        xfrac = round(x_o + strength * dx, 4)

    return xfrac

def mutate(run_id, generation):
    print( "\nCreating generation :\t%s" % (generation) )
    wd = os.environ['HTSOHM_DIR']
    md = os.environ['MAT_DIR']
    fd = os.environ['FF_DIR']
    # Load parameter boundaries from ../config/<run_data>.yaml
    with open(wd + '/config/' + run_id + '.yaml') as yaml_file:
        config = yaml.load(yaml_file)
    
    children_per_generation = config["children-per-generation"]
    number_of_atomtypes = config["number-of-atom-types"]
    
    number_density_limits = config["number-density-limits"]
    ndenmin = number_density_limits[0]
    ndenmax = number_density_limits[1]

    lattice_constant_limits = config["lattice-constant-limits"]
    xmin = lattice_constant_limits[0]
    xmax = lattice_constant_limits[1]
    ymin = lattice_constant_limits[0]
    ymax = lattice_constant_limits[1]
    zmin = lattice_constant_limits[0]
    zmax = lattice_constant_limits[1]

    epsilon_limits = config["epsilon-limits"]
    epmin = epsilon_limits[0]
    epmax = epsilon_limits[1]

    sigma_limits = config["sigma-limits"]
    sigmin = sigma_limits[0]
    sigmax = sigma_limits[1]

    eq = config["elemental-charge"]

    first = generation * children_per_generation
    last = (generation + 1) * children_per_generation
    child_ids = np.arange(first, last)
    strength_array = np.load(wd + '/config/' + run_id + '.npy')   # Load strength-parameter array

    for i in child_ids:
        child_id = str(i)
        child = session.query(RunData).filter(
            RunData.run_id == run_id, RunData.material_id == child_id)
        for item in child:
            parent_id = item.parent_id
        parent = session.query(RunData).get(parent_id)
        parent_material_id = parent.material_id
        parent_ml_bin = parent.methane_loading_bin
        parent_sa_bin = parent.surface_area_bin
        parent_vf_bin = parent.void_fraction_bin
        strength = strength_array[parent_ml_bin, parent_sa_bin, parent_vf_bin]
        
        pd = "%s/%s-%s" % (fd, run_id, parent_material_id) # Parent's forcefield directory
        cd = "%s/%s-%s" % (fd, run_id, child_id)           # Child's forcefield directory
        os.mkdir(cd)

        # Copy force_field.def
        shutil.copy("%s/force_field.def" % (pd), cd)
        shutil.copy("%s/pseudo_atoms.def" % (pd), cd)      # CHANGE THIS WHEN ADDING CHARGES

        # Load data from parent's definition files
        parent_mixing = "%s/force_field_mixing_rules.def" % (pd)
        n1, n2, ep_o, sig_o = np.genfromtxt(parent_mixing, unpack=True,
            skip_header=7, skip_footer=9)
        parent_cif = "%s/%s-%s.cif" % (md, run_id, parent_material_id)
        cif_atype = np.genfromtxt(parent_cif, usecols=0, dtype=str, 
            skip_header=16)
        n1, n2, x_o, y_o, z_o = np.genfromtxt(parent_cif, unpack=True,
            skip_header=16)
        a_foot = len(x_o) + 11
        b_foot = len(x_o) + 10
        c_foot = len(x_o) + 9
        n1, a_o = np.genfromtxt(parent_cif, unpack=True, skip_header=4,
            skip_footer=a_foot)
        n1, b_o = np.genfromtxt(parent_cif, unpack=True, skip_header=5,
            skip_footer=b_foot)
        n1, c_o = np.genfromtxt(parent_cif, unpack=True, skip_header=6,
            skip_footer=c_foot)

        # Open child definition files
        cif_file = open("%s/%s-%s.cif" % (md, run_id, child_id), "w")
        mix_file = open("%s/force_field_mixing_rules.def" % (cd), "w")

        # Perturb crystal lattice parameters, then write to file
        a_r = round(random() * (xmax - xmin) + xmin, 4)
        b_r = round(random() * (ymax - ymin) + ymin, 4)
        c_r = round(random() * (zmax - zmin) + zmin, 4)
        a_n = round(a_o + strength * (a_r - a_o), 4)
        b_n = round(b_o + strength * (b_r - b_o), 4)
        c_n = round(c_o + strength * (c_r - c_o), 4)

        cif_file.write(
            "\nloop_\n" +
            "_symmetry_equiv_pos_as_xyz\n" +
            "  x,y,z\n" +
            "_cell_length_a\t%s\n" % (a_n) +
            "_cell_length_b\t%s\n" % (b_n) +
            "_cell_length_c\t%s\n" % (c_n) +
            "_cell_angle_alpha\t90.0000\n" +
            "_cell_angle_beta\t90.0000\n" +
            "_cell_angle_gamma\t90.0000\n" +
            "loop_\n" +
            "_atom_site_label\n" +
            "_atom_site_type_symbol\n" +
            "_atom_site_fract_x\n" +
            "_atom_site_fract_y\n" +
            "_atom_site_fract_z\n")

        # Perturb number density
        n_o = len(x_o)
        nden_o = n_o / (a_o * b_o * c_o)
        nden_r = round(random() * (ndenmax - ndenmin) + ndenmin, 4)
        nden = nden_o + strength * (nden_r - nden_o)
        n = int(nden * a_n * b_n * c_n)

        # Remove pseudo-atoms from unit cell (if necessary)...
        omit_n = 0
        if n < n_o:
            omit_n = n_o - n
        ndiff = n_o - n  # rewrite this in a more logical

        # Perturbing atomic positions (xfrac, yfrac, zfrac)...
        for l in range(n_o - omit_n):
            x_r = random()
            y_r = random()
            z_r = random()
            x = delta_x(x_o[l], x_r, strength)
            y = delta_x(y_o[l], y_r, strength)
            z = delta_x(z_o[l], z_r, strength)
            charge = 0.  # if no charge in .cif---ERRRORRR
            cif_file.write("%s\tC\t%s\t%s\t%s\n" % (cif_atype[l], x, y, z)
        # Add pseudo-atoms to unit cell (if necessary)...
        if n > n_o:
            add_n = n - n_o
            for m in range(add_n):
                atyp = choice(cif_atype)
                charge = 0.  # CHANGE THIS!!!
                x = round(random(), 4)
                y = round(random(), 4)
                z = round(random(), 4)
                cif_file.write( "%s\tC\t%s\t%s\t%s\n" % (cif_atype[l], x, y, z)

        mix_file.write(
            "# general rule for shifted vs truncated\n" +
            "shifted\n" +
            "# general rule tailcorrections\n" +
            "no\n" +
            "# number of defined interactions\n" +
            "%s\n" % (number_of_atomtypes + 8) +
            "# type interaction, parameters.    IMPORTANT:" +
            " define shortest matches first, so that more " +
            "specific ones overwrites these\n")

        # Perturbing LJ parameters (sigma, epsilon)...
        for o in range(number_of_atomtypes):
            ep_r = round(random() * (epmax - epmin) + epmin, 4)
            sig_r = round(random() * (sigmax - sigmin) + sigmin, 4)
            epsilon = round(ep_o[o] + strength * (ep_r - ep_o[o]), 4)
            sigma = round(sig_o[o] + strength * (sig_r - sig_o[o]), 4)
            row = "A_%s\tlennard-jones\t%s\t%s\n" % (o, epsilon, sigma))
            mix_file.write(row)

        mix_file.write(
            "N_n2\tlennard-jones\t36.0\t3.31\n" +
            "N_com\tnone\n" +
            "C_co2\tlennard-jones\t27.0\t2.80\n" +
            "O_co2\tlennard-jones\t79.0\t3.05\n" +
            "CH4_sp3\tlennard-jones\t158.5\t3.72\n" +
            "He\tlennard-jones\t10.9\t2.64\n" +
            "H_h2\tnone\n" +
            "H_com\tlennard-jones\t36.7\t2.958\n" +
            "# general mixing rule for Lennard-Jones\n" +
            "Lorentz-Berthelot")
        cif_file.close()
        mix_file.close()

    print( "...done!\n" )
