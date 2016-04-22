import os
import shutil

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
import numpy as np
from random import random, choice

from runDB_declarative import Base, RunData
import binning as bng
from simulate import GetValue, id_to_mat


# Create "strength" array...
def FirstS(run_ID, strength_0):

    bins = int(run_ID[-1])

    s_array = np.zeros([bins, bins, bins])

    for i in range(bins):
        for j in range(bins):
            for k in range(bins):
                s_array[i, j, k] = strength_0

    wd = os.environ['HTSOHM_DIR']
    np.save(wd + '/' + run_ID, s_array)


def CalculateS(run_ID, generation):

    wd = os.environ['HTSOHM_DIR']

    with open(wd + '/' + run_ID + '.txt') as origin:
        for line in origin:
            if "Children per generation:" in line:
                children_per_generation = int(line.split()[3])
    First = generation * children_per_generation
    Last = (generation + 1) * children_per_generation
    c_IDs = np.arange(First, Last)
    gen0_IDs = np.arange(0, children_per_generation)

    Strength_o = np.load(wd + '/' + run_ID + '.npy')       # Load strength-parameter array

    p_IDs = []
    p_bins = []
    for i in c_IDs:
        p_ID = GetValue(run_ID, i, "Mat")
        p_IDs.append( p_ID )
        p_MLb = GetValue(run_ID, p_ID, "Bin_ML")
        p_SAb = GetValue(run_ID, p_ID, "Bin_SA")
        p_VFb = GetValue(run_ID, p_ID, "Bin_VF")

        p_bin = [ p_MLb, p_SAb, p_VFb ]
        p_bins.append( p_bin )

    pgen_IDs = np.arange( First - children_per_generation,
                          Last - children_per_generation )
    gp_IDs = []
    gp_bins = []
    for i in pgen_IDs:
        gp_ID = GetValue(run_ID, i, "id")
        gp_IDs.append( gp_ID )
        gp_MLb = GetValue(run_ID, gp_ID, "Bin_ML")
        gp_SAb = GetValue(run_ID, gp_ID, "Bin_SA")
        gp_VFb = GetValue(run_ID, gp_ID, "Bin_VF")
        gp_bin = [ gp_MLb, gp_SAb, gp_VFb ]
        gp_bins.append( gp_bin )
    
    bin_list = []
    for i in p_bins:
        if i in gp_bins:
            bin_list.append(i)

    dS_bins = []
    for i in bin_list:
        if i not in dS_bins:
            dS_bins.append(i)

    counts = bng.CountAll(run_ID)

    bin_counts = []
    for i in dS_bins:
        p_bin = i
        p_count = counts[ i[0],i[1],i[2] ]
        
        c_bins = []
        c_counts = []
        for j in range(len(gp_bins)):
            if p_bin == gp_bins[j]:
                ID = j + children_per_generation
                ML_bin = GetValue(run_ID, ID, "Bin_ML")
                SA_bin = GetValue(run_ID, ID, "Bin_SA")
                VF_bin = GetValue(run_ID, ID, "Bin_VF")
                c_bin = [ ML_bin, SA_bin, VF_bin ]

                if c_bin not in c_bins:
                    c_bins.append(c_bin)

                    count = int( counts[ c_bin[0], c_bin[1], c_bin[2] ] )
                    c_counts.append(count)

        p_data = [p_bin, p_count]
        c_data = [c_bins, c_counts]
        row = [p_data, c_data]
        bin_counts.append(row)

    for i in bin_counts:

        p_bin = i[0][0]
        p_count = i[0][1]

        c_bins = i[1][0]
        c_counts = i[1][1]

        a = p_bin[0]
        b = p_bin[1]
        c = p_bin[2]

        S_0 = Strength_o[a, b, c]

        if p_bin not in c_bins:
            Strength_o[a, b, c] = 0.5 * S_0

        if p_bin in c_bins:

            if p_count < 1.1 * min(c_counts):
                Strength_o[a, b, c] = 0.5 * S_0

            pos = c_bins.index(p_bin)
            val = 0
            for j in range(len(c_bins)):
                if j != pos:
                    if c_counts[j] > val:
                        val = c_counts[j]

            if p_count >= 3 * val:
                Strength_o[a, b, c] = 1.5 * S_0

    S_file = "%s/%s.npy" % (wd, run_ID)
    os.remove(S_file)
    np.save(S_file, Strength_o)


# function for finding "closest" distance over periodic boundaries
def closestDist(x_o, x_r):
    a = 1 - x_r + x_o
    b = abs(x_r - x_o)
    c = 1 - x_o + x_r
    dx = min(a, b, c)

    return dx # is closest in R1 the closest in R3???


# given intial and "random" x-fraction, returns new x fraction
def deltax(x_o, x_r, strength):  # removed random()
    dx = closestDist(x_o, x_r)

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


def mutate(run_ID, generation):

    print( "\nCreating generation :\t%s" % (generation) )

    wd = os.environ['HTSOHM_DIR']
    md = os.environ['MAT_DIR']
    fd = os.environ['FF_DIR']

    # Load parameter boundaries from `run_ID`.txt
    with open(wd + '/' + run_ID + '.txt') as origin:
        for line in origin:
            if "Children per generation:" in line:
                children_per_generation = int(line.split()[3])
            elif "Number of atom-types:" in line:
                number_of_atomtypes = int(line.split()[3])
            elif "Number density:" in line:
                ndenmin = float(line.split()[2])
                ndenmax = float(line.split()[4])
            elif "Lattice constant, a:" in line:
                xmin = float(line.split()[3])
                xmax = float(line.split()[5])
            elif "Lattice constant, b:" in line:
                ymin = float(line.split()[3])
                ymax = float(line.split()[5])
            elif "Lattice constant, c:" in line:
                zmin = float(line.split()[3])
                zmax = float(line.split()[5])
            elif "Epsilon:" in line:
                epmin = float(line.split()[1])
                epmax = float(line.split()[3])
            elif "Sigma:" in line:
                sigmin = float(line.split()[1])
                sigmax = float(line.split()[3])
            elif "Elemental charge:" in line:
                eq = float(line.split()[2])

    first = generation * children_per_generation
    last = (generation + 1) * children_per_generation
    child_IDs = np.arange(first, last)
    
    Strength = np.load(wd + '/' + run_ID + '.npy')         # Load strength-parameter array


    for i in child_IDs:
        child_ID = str(i)
        p = GetValue(run_ID, child_ID, "Parent")        # Find parent ID
        p_ID = id_to_mat(run_ID, p)
        p_MLb = GetValue(run_ID, p_ID, "Bin_ML")           # Find parent-bin coordinates
        p_SAb = GetValue(run_ID, p_ID, "Bin_SA")
        p_VFb = GetValue(run_ID, p_ID, "Bin_VF")
        strength = Strength[p_MLb, p_SAb, p_VFb]
        
        pd = "%s/%s-%s" % (fd, run_ID, p_ID)               # Parent's forcefield directory
        cd = "%s/%s-%s" % (fd, run_ID, child_ID)           # Child's forcefield directory
        os.mkdir(cd)

        # Copy force_field.def
        shutil.copy("%s/force_field.def" % (pd), cd)
        shutil.copy("%s/pseudo_atoms.def" % (pd), cd) # CHANGE THIS WHEN ADDING CHARGES

        # Load data from parent's definition files
        n1, n2, ep_o, sig_o = np.genfromtxt("%s/force_field_mixing_rules.def" % (pd),
                                            unpack=True,
                                            skip_header=7,
                                            skip_footer=9)
        cif_atype = np.genfromtxt("%s/%s-%s.cif" % (md, run_ID, p_ID),
                                  usecols=0, dtype=str, skip_header=16)
        n1, n2, x_o, y_o, z_o = np.genfromtxt("%s/%s-%s.cif" % (md, run_ID, p_ID),
                                              unpack=True, skip_header=16)
        a_foot = len(x_o) + 11
        b_foot = len(x_o) + 10
        c_foot = len(x_o) + 9
        n1, a_o = np.genfromtxt("%s/%s-%s.cif" % (md, run_ID, p_ID),
                                unpack=True, skip_header=4,
                                skip_footer=a_foot)
        n1, b_o = np.genfromtxt("%s/%s-%s.cif" % (md, run_ID, p_ID),
                                unpack=True, skip_header=5,
                                skip_footer=b_foot)
        n1, c_o = np.genfromtxt("%s/%s-%s.cif" % (md, run_ID, p_ID),
                                unpack=True, skip_header=6,
                                skip_footer=c_foot)

        # Open child definition files
        cif_file = open("%s/%s-%s.cif" % (md, run_ID, child_ID), "w")
        mix_file = open("%s/force_field_mixing_rules.def" % (cd), "w")

        # Perturb crystal lattice parameters, then write to file
        a_r = round(random() * (xmax - xmin) + xmin, 4)
        b_r = round(random() * (ymax - ymin) + ymin, 4)
        c_r = round(random() * (zmax - zmin) + zmin, 4)
        a_n = round(a_o + strength * (a_r - a_o), 4)
        b_n = round(b_o + strength * (b_r - b_o), 4)
        c_n = round(c_o + strength * (c_r - c_o), 4)

        cif_file.write("\nloop_\n" +
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

            xfrac = deltax(x_o[l], x_r, strength)
            yfrac = deltax(y_o[l], y_r, strength)
            zfrac = deltax(z_o[l], z_r, strength)

            charge = 0.  # if no charge in .cif---ERRRORRR
            
            cif_file.write( "%s\tC\t%s\t%s\t%s\n" % (cif_atype[l], xfrac, yfrac, zfrac) )

        # Add pseudo-atoms to unit cell (if necessary)...
        if n > n_o:
            add_n = n - n_o
            for m in range(add_n):

                atyp = choice(cif_atype)

                charge = 0.  # CHANGE THIS!!!

                xfrac = round(random(), 4)
                yfrac = round(random(), 4)
                zfrac = round(random(), 4)

                cif_file.write( "%s\tC\t%s\t%s\t%s\n" % (cif_atype[l], xfrac, yfrac, zfrac) )

        mix_file.write( "# general rule for shifted vs truncated\n" +
                        "shifted\n" +
                        "# general rule tailcorrections\n" +
                        "no\n" +
                        "# number of defined interactions\n" +
                        "%s\n" % (number_of_atomtypes + 8) +
                        "# type interaction, parameters.    IMPORTANT:" +
                          " define shortest matches first, so that more " +
                          "specific ones overwrites these\n" )

        # Perturbing LJ parameters (sigma, epsilon)...
        for o in range(number_of_atomtypes):

            ep_r = round(random() * (epmax - epmin) + epmin, 4)
            sig_r = round(random() * (sigmax - sigmin) + sigmin, 4)

            epsilon = round(ep_o[o] + strength * (ep_r - ep_o[o]), 4)
            sigma = round(sig_o[o] + strength * (sig_r - sig_o[o]), 4)

            mix_file.write( "A_%s\tlennard-jones\t%s\t%s\n" % (o, epsilon, sigma) )

        mix_file.write( "N_n2\tlennard-jones\t36.0\t3.31\n" +
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
