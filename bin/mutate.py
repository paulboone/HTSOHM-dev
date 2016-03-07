import os
import numpy as np
from binning import *
import shutil
from random import random, choice

def firstS(p_dir, m_strength, bins):

    s_array = np.zeros( [bins, bins, bins] )

    for i in range(bins):
        for j in range(bins):
            for k in range(bins):
                s_array[i, j , k] = m_strength

    np.save(p_dir + '/s_file', s_array)
    

def calc_S(p_dir, bins): 

    gp_dir = 'gen' + str(int(p_dir[-1]) - 1)

    gp_S = np.load( gp_dir + '/s_file.npy' )

    p_list = np.genfromtxt(p_dir + '/p_list.txt', usecols=0, dtype=str)

    p_Xs = np.genfromtxt(p_dir + '/p_list.txt', usecols=1, dtype=str)
    p_Ys = np.genfromtxt(p_dir + '/p_list.txt', usecols=2, dtype=str)
    p_Zs = np.genfromtxt(p_dir + '/p_list.txt', usecols=3, dtype=str)

    p_bins = []
    for i in range( len(p_Xs) ):
        pos = [ int( p_Xs[i][1:-1] ),
                int( p_Ys[i][:-1] ),
                int( p_Zs[i][:-1] ) ]
        p_bins.append(pos)

    gp_list = np.genfromtxt(gp_dir + '/p_list.txt', usecols=0, dtype=str)

    gp_Xs = np.genfromtxt(gp_dir + '/p_list.txt', usecols=1, dtype=str)
    gp_Ys = np.genfromtxt(gp_dir + '/p_list.txt', usecols=2, dtype=str)
    gp_Zs = np.genfromtxt(gp_dir + '/p_list.txt', usecols=3, dtype=str)

    gp_bins = []
    for i in range( len(gp_Xs) ):
        pos = [ int( gp_Xs[i][1:-1] ),
                int( gp_Ys[i][:-1] ),
                int( gp_Zs[i][:-1] ) ]
        gp_bins.append(pos)

    bin_list = []
    for i in p_bins:
        if i in gp_bins:
            bin_list.append(i)
        
    dS_bins = []
    for i in bin_list:
        if i not in dS_bins:
            dS_bins.append(i)
        
    number = len( np.genfromtxt(gp_dir + '/ch4_abs_cc_cc.txt',
                                usecols=0, dtype=str) )

    counts, ID_array = bin3d(p_dir[1:], bins)

    bin_counts = []

    for i in dS_bins:
    
        p_bin = i
    
        p_count = counts[ i[0],
                          i[1],
                          i[2] ]
    
        c_bins = []
        c_counts = []
    
        for j in range( len(gp_bins) ):
        
            if p_bin == gp_bins[j]:
            
                ID = j + number
                child = 'MAT-' + str(ID)
            
                c_bin = find_bin(child, ID_array)
            
                if c_bin not in c_bins:
                    c_bins.append(c_bin)
                
                    count = int(counts[ bin_c[0],
                                        bin_c[1],
                                        bin_c[2] ])
                
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
    
        S_0 = gp_S[a, b, c]
    
        if p_bin not in c_bins:
            gp_S[a, b, c] = 0.5 * S_0
        
        if p_bin in c_bins:
        
            if p_count < 1.1 * min(c_counts):
                gp_S[a, b, c] = 0.5 * S_0
        
            pos = c_bins.index(p_bin)
            val = 0
            for j in range( len(c_bins)):
                if j != pos:
                    if c_counts[j] > val:
                        val = c_counts[j]
        
            if p_count >= 3 * val:
                gp_S[a, b, c] = 1.5 * S_0

    np.save(p_dir + '/s_file', gp_S)

# function for finding "closest" distance over periodic boundaries
def closestDist(x_o, x_r):  
    a = 1 - x_r + x_o
    b = abs(x_r - x_o)
    c = 1 - x_o + x_r
    dx = min(a, b, c)

    return dx
#is closest in R1 the closest in R3???

# given intial and "random" x-fraction, returns new x fraction
def deltax(x_o, x_r, strength): # removed random()
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

def mutate(p_dir,     #directory containing parent library
           n_atype,   #number of atom-types
           c_dir):     #directory for new, child library

    #boundaries
    xmin = 25.6
    xmax = 51.2
    ymin = xmin
    ymax = xmax
    zmin = xmin
    zmax = xmax
    ndenmin = 0.0000013905
    ndenmax = 0.04302
    epmin = 1.2580
    epmax = 513.264
    sigmin = 1.052342
    sigmax = 6.549291

    p_list = np.genfromtxt(p_dir + '/p_list.txt', usecols=0, dtype=str)
    p_Xs = np.genfromtxt(p_dir + '/p_list.txt', usecols=1, dtype=str)
    p_Ys = np.genfromtxt(p_dir + '/p_list.txt', usecols=2, dtype=str)
    p_Zs = np.genfromtxt(p_dir + '/p_list.txt', usecols=3, dtype=str)
    
    p_bins = []
    for i in range( len(p_Xs) ):
        pos = [ int( p_Xs[i][1:-1] ),
                int( p_Ys[i][:-1] ),
                int( p_Zs[i][:-1] ) ]
        p_bins.append( pos )

    s_array = np.load( p_dir + '/s_file.npy' )

#    s_list = np.genfromtxt(p_dir + '/s_list.txt')

    top_dir = c_dir

    if not os.path.exists(top_dir):
        os.mkdir(top_dir)

    n_child = len(p_list)
    print( 'making ' + str(n_child) + ' materials...' )

    p_lib = np.genfromtxt(p_dir + '/ch4_abs_cc_cc.txt', usecols=0, dtype=str)
    last_mat = int(p_lib[-1][4:])

    for i in np.arange(last_mat + 1, last_mat + n_child + 1):

        c_ID = 'MAT-' + str(i)

        p_num = i - last_mat - 1
        p_ID = p_list[p_num]

        loc = p_bins[p_num]

        strength = s_array[ loc[0], loc[1], loc[2] ]
  
        # creating child directory...
        mat_dir = top_dir + '/' + c_ID
        os.mkdir(mat_dir)

        # copying parent data files...
        shutil.copy(p_dir + '/' + p_ID + 
                    '/pseudo_atoms.def', mat_dir)  # consider charges!
        shutil.copy(p_dir + '/' + p_ID + 
                    '/force_field.def', mat_dir)

        # importing values from parent data files...
        n1, n2, ep_o, sig_o = np.genfromtxt(p_dir + '/' + p_ID + 
                                            '/force_field_mixing_rules.def',
                                            unpack=True,
                                            skip_header=7,
                                            skip_footer=9)

        cif_atyp = np.genfromtxt(p_dir + '/' + p_ID  + '/' + p_ID + '.cif',
                                 usecols=0, dtype=str, skip_header=17)

        n1, n2, x_o, y_o, z_o = np.genfromtxt(p_dir + '/' + p_ID + 
                                              '/' + p_ID + '.cif',
                                              unpack=True, skip_header=17)

        # opening child data files ...            
        cif_file = open(os.path.abspath(mat_dir) + '/' + c_ID + '.cif', 'w')

        plog_old = p_dir + '/' + p_ID + '/parents_log.txt'

        if not os.path.isfile(plog_old):
            parents_log = open(os.path.abspath(mat_dir) + '/parents_log.txt', 'w')
            parents_log.write(p_ID + '\n' + c_ID)
        if os.path.isfile(plog_old):
            plog = np.genfromtxt(plog_old, dtype=str)
            parents_log = open(os.path.abspath(mat_dir) + '/parents_log.txt', 'w')
            for f in plog:
                parents_log.write(f + '\n')
            parents_log.write(c_ID)

        mixing_file = open(os.path.abspath(mat_dir) + 
                           '/force_field_mixing_rules.def', 'w')

        # perturbing crystal lattice parameters (a,b,c)...        
        a_footer = len(x_o)+11
        b_footer = len(x_o)+10
        c_footer = len(x_o)+9

        n1, a_o = np.genfromtxt(p_dir + '/' + p_ID + '/' + p_ID + '.cif', 
                                unpack=True, skip_header=5, skip_footer=a_footer)

        n1, b_o = np.genfromtxt(p_dir + '/' + p_ID + '/' + p_ID  + '.cif',
                                unpack=True, skip_header=6, skip_footer=b_footer)

        n1, c_o = np.genfromtxt(p_dir + '/' + p_ID + '/' + p_ID + '.cif',
                                unpack=True, skip_header=7, skip_footer=c_footer)

        a_r = round(random() * (xmax - xmin) + xmin, 4)
        b_r = round(random() * (ymax - ymin) + ymin, 4)
        c_r = round(random() * (zmax - zmin) + zmin, 4)

        a_n = round(a_o + strength * (a_r - a_o), 4)
        b_n = round(b_o + strength * (b_r - b_o), 4)
        c_n = round(c_o + strength * (c_r - c_o), 4)

        # writing new crystal lattice parameters to file...                    
        cif_head = (c_ID +
                    '\n\nloop_\n' +
                    '_symmetry_equiv_pos_as_xyz\n' +
                    '  x,y,z\n' +
                    '_cell_length_a          ' + str(a_n) + '\n' +
                    '_cell_length_b          ' + str(b_n) + '\n' +
                    '_cell_length_c          ' + str(c_n) + '\n' +
                    '_cell_angle_alpha       90.0000\n' +
                    '_cell_angle_beta        90.0000\n' +
                    '_cell_angle_gamma       90.0000\n' +
                    'loop_\n' +
                    '_atom_site_label\n' +
                    '_atom_site_type_symbol\n' +
                    '_atom_site_fract_x\n' +
                    '_atom_site_fract_y\n' +
                    '_atom_site_fract_z\n')

        cif_file.write(cif_head)

        # perturbing number density...
        n_o = len(x_o)
        nden_o = n_o / (a_o * b_o * c_o)

        nden_r = round(random() * (ndenmax - ndenmin) + ndenmin, 4)

        nden = nden_o + strength * (nden_r - nden_o)
        n = int(nden * a_n * b_n * c_n)

        # removing pseudo-atoms from unit cell (if necessary)...
        omit_n = 0
        if n < n_o:
            omit_n = n_o - n

        ndiff = n_o - n # rewrite this in a more logical

        # perturbing atomic positions (xfrac, yfrac, zfrac)...
        for l in range(n_o - omit_n):

            x_r = random()
            y_r = random()
            z_r = random()

            xfrac = deltax(x_o[l], x_r, strength)
            yfrac = deltax(y_o[l], y_r, strength)
            zfrac = deltax(z_o[l], z_r, strength)

            charge = 0. #if no charge in .cif---ERRRORRR

            cif_line = (cif_atyp[l] + '     C     ' + str(xfrac) +
                        '     ' + str(yfrac) + '     ' + 
                        str(zfrac) + '\n')
            cif_file.write(cif_line)

        # adding pseudo-atoms to unit cell (if necessary)...
        if n > n_o:
            add_n = n - n_o
            for m in range(add_n):

                atyp = choice(cif_atyp)

                charge = 0. # CHANGE THIS!!!                           

                xfrac = round(random(), 4)
                yfrac = round(random(), 4)
                zfrac = round(random(), 4)

                new_line = (cif_atyp[l] + '     C     ' + str(xfrac) +
                            '     ' + str(yfrac) + '     ' + 
                            str(zfrac) + '\n')
                cif_file.write(new_line)

        mixing_head = ('# general rule for shifted vs truncated\n' +
                       'shifted\n' +
                       '# general rule tailcorrections\n' +
                       'no\n' +
                       '# number of defined interactions\n' +
                       str(n_atype + 8) + '\n' +
                       '# type interaction, parameters.    IMPORTANT:' + 
                       ' define shortest matches first, so that more ' +
                       'specific ones overwrites these\n')
        mixing_file.write(mixing_head)

        # perturbing LJ parameters (sigma, epsilon)...
        for o in range(n_atype):

            ep_r = round(random() * (epmax - epmin) + epmin, 4)
            sig_r = round(random() * (sigmax - sigmin) + sigmin, 4)

            epsilon = round(ep_o[o] + strength * (ep_r - ep_o[o]), 4) 
            sigma = round(sig_o[o] + strength * (sig_r - sig_o[o]), 4)

            mixing_line = ('A_' + str(o) + '   lennard-jones   ' + 
                                       str(epsilon) + '   ' + str(sigma) + '\n')
            mixing_file.write(mixing_line)

        mixing_foot = ('N_n2        lennard-jones   36.0     3.31\n' +
                       'N_com       none\n' +
                       'C_co2       lennard-jones   27.0     2.80\n' +
                       'O_co2       lennard-jones   79.0     3.05\n' +
                       'CH4_sp3     lennard-jones   158.5    3.72\n' +
                       'He          lennard-jones   10.9     2.64\n' +
                       'H_h2        none\n' +
                       'H_com       lennard-jones   36.7     2.958\n' +
                       '# general mixing rule for Lennard-Jones\n' +
                       'Lorentz-Berthelot')
        mixing_file.write(mixing_foot)

        cif_file.close() 
        mixing_file.close()
    print( '...done!' )
