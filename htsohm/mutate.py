import os
import shutil
import yaml

import numpy as np
from random import random, choice

from htsohm import binning as bng
from htsohm import simulate as sim
from htsohm.runDB_declarative import session, RunData

def create_strength_array(run_id):
    """Create 3-dimensional array of strength parameters.
    The method divides the 3-dimensional materials parameter space into congruent bins. Each bin
    has its own mutation strength, which is automatically increased or decreased depending upon
    the distribution of children whose parent belong to that bin. This function initializes a
    3-dimensional strength-parameter bin with the initial mutation strength chosen for the run.
    """
    wd = os.environ['HTSOHM_DIR']
    config_file = os.path.join(wd, 'config', run_id + '.yaml')
    with open(config_file) as yaml_file:
        config = yaml.load(yaml_file)
    bins               = config["number-of-bins"]
    initial_strength   = config["initial-mutation-strength"]
    strength_array = initial_strength * np.ones([bins, bins, bins])
    filename = os.path.join(wd, 'config', run_id)
    np.save(strength_array, filename)

def recalculate_strength_array(run_id, generation):
    """Rewrite 3-dimensional array of stregnth parameters.
    In the event that a particular bin contains parents whose children exhibit radically
    divergent properties, the strength paramter for the bin is modified. In order to deftermine
    which bins to adjust, the script refers to the distribution of children in the previous
    generation which share a common parent. The criteria follows:
     ________________________________________________________________
     - if none of the children share  |  decrease strength parameter
       the parent's bin               |  by 50 %
     - if all other child-bin-counts  |  
       are at least 10% greater than  |
       the parent-bin-count           |
     _________________________________|_____________________________
     - if the parent-bin-count is at  |  increase strength parameter
       least 3 times greater than     |  by 50 %
       all other child-bin-counts     |
     _________________________________|_____________________________
    """
    s = session
    db = RunData
    
    parent_list = []
    children = s.query(db).filter(db.run_id == run_id, db.generation == generation)
    for item in children:
        p = s.query(db).get(item.parent_id)
        parent = {
            "id"  : item.parent_id,
            "bin" : [p.methane_loading_bin, p.surface_area_bin, p.void_fraction_bin]}
        parent_list.append(parent)

    wd = os.environ['HTSOHM_DIR']
    strength_array_file = os.path.join(wd, 'config', run_id + '.npy')
    strength_array = np.load(strength_array_file)

    for parent in parent_list:
        parent_id    = parent["id"]
        parent_bin   = parent["bin"]
        child_bins   = []
        child_counts = []
        children = s.query(db).filter(db.run_id == run_id, db.generation == generation - 1,
            db.parent_id == parent_id).all()
        for child in children:
            child_bin = [child.methane_loading_bin, child.surface_area_bin, child.void_fraction_bin]
            bin_count = s.query(db).filter(db.run_id == run_id,
                db.generation == generation - 1, db.parent_id == parent_id,
                db.methane_load_bin == child[0], db.surface_area_bin == child[1],
                db.void_fraction_bin == child[2]).count()
            if child_bin not in child_bins:
                child_bins.append(child_bin)
                child.counts.append(bin_count)
        if parent_bin not in child_bins:
            strength_array[parent_bin] = 0.5 * strength_array[parent_bin]
        elif parent_bin in child_bins:
            parent_bin_count = child_counts[child_bins.index(parent_bin)]
            del child_counts[child_bins.index(parent_bin)]
            if parent_bin_count < 1.1 * min(child_counts):
                strength_array[parent_bin] = 0.5 * strength_array[parent_bin]
            if parent_bin_count >= 3 * min(child_counts):
                strength_array[parent_bin] = 1.5 * strength_array[parent_bin]

    os.remove(strength_array_file)                   # remove old file
    np.save(strength_array_file, strength_array)     # write new file

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
    config_file = os.path.join(wd, 'config', run_id + '.yaml')
    with open(config_file) as yaml_file:
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
        
        pd = os.path.join(fd, run_id + '-' + str(parent_material_id)) # parent FF directory
        cd = os.path.join(fd, run_id + '-' + str(child_id))           # child FF directory 
        os.mkdir(cd)

        # Copy force_field.def
        parent_force = os.path.join(pd, 'force_field.def')
        shutil.copy(parent_force, cd)
        parent_pseudo = os.path.join(pd, 'pseudo_atoms.def')
        shutil.copy(parent_pseudo, cd)

        # Load data from parent's definition files
        parent_mixing = os.path.join(pd, 'force_field_mixing_rules.def')
        n1, n2, ep_o, sig_o = np.genfromtxt(parent_mixing, unpack=True,
            skip_header=7, skip_footer=9)
        parent_cif_filename = "%s-%s.cif" % (run_id, parent_material_id)
        parent_cif = os.path.join(md, parent_cif_filename)
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
        child_cif_filename = "%s-%s.cif" % (run_id, child_id)
        cif_file = os.path.join(md, child_cif_filename)
        mix_file = os.path.join(cd, 'force_field_mixing_rules.def')

        # Perturb crystal lattice parameters, then write to file
        a_r = round(random() * (xmax - xmin) + xmin, 4)
        b_r = round(random() * (ymax - ymin) + ymin, 4)
        c_r = round(random() * (zmax - zmin) + zmin, 4)
        a_n = round(a_o + strength * (a_r - a_o), 4)
        b_n = round(b_o + strength * (b_r - b_o), 4)
        c_n = round(c_o + strength * (c_r - c_o), 4)

        with open(cif_file, "w") as file:
            file.write(
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
            with open(cif_file, "a") as file:
                file.write("%s\tC\t%s\t%s\t%s\n" % (cif_atype[l], x, y, z))
        # Add pseudo-atoms to unit cell (if necessary)...
        if n > n_o:
            add_n = n - n_o
            for m in range(add_n):
                atyp = choice(cif_atype)
                charge = 0.  # CHANGE THIS!!!
                x = round(random(), 4)
                y = round(random(), 4)
                z = round(random(), 4)
                with open(cif_file, "a") as file:
                    file.write( "%s\tC\t%s\t%s\t%s\n" % (cif_atype[l], x, y, z))

        with open(mix_file, "w") as file:
            file.write(
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
            row = "A_%s\tlennard-jones\t%s\t%s\n" % (o, epsilon, sigma)
            with open(mix_file, "a") as file:
                file.write(row)

        with open(mix_file, "a") as file:
            file.write(
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

    print( "...done!\n" )
