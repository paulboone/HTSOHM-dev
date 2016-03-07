import numpy as np
import pylab as py
import os
from numpy.random import choice #python 2.7.8... not earlier or py3
from random import randrange, random, choice
import shutil

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

def mutate(data_path,
           atom_types,
           number_of_children,
           strength,
           bins,
           bin_dimensions=['HV', 'SA', 'CH4']):
    
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

    if 'HV' in bin_dimensions:

        HV_labels = np.genfromtxt(data_path + '/HVdata2col.txt', usecols=0,
                                      dtype=str)
        HV_values = np.genfromtxt(data_path + '/HVdata2col.txt', usecols=1,
                                      dtype=float)

        HV_binmin = 0.
        HV_binmax = 1.
        HV_step = HV_binmax/bins
        edges_HV = np.arange(HV_binmin, HV_binmax + HV_step, HV_step)

        HV_IDs = np.empty(bins).tolist()
        for i in range(bins):
            HV_IDs[i] = []

        # 4 debug
        HV_vals = np.empty(bins).tolist()
        for i in range(bins):
            HV_vals[i] = []

        counter=0
        for i in range(bins):
            for j in range(len(HV_labels)):
                if HV_values[j] >= edges_HV[i] and HV_values[j] < edges_HV[i + 1]:
                    counter+=1
                    HVID = HV_labels[j]
                    HV_IDs[i].append(HVID)
                    HV_vals[i].append(HV_values[j])

    if 'SA' in bin_dimensions:

        SA_labels = np.genfromtxt(data_path + '/SAdata_m2_cc.txt', usecols=0,
                                  dtype=str)
        SA_values = np.genfromtxt(data_path + '/SAdata_m2_cc.txt', usecols=1,
                                  dtype=float)

        SA_binmin = 0.
        SA_binmax = 4500.
        SA_step = SA_binmax/bins
        edges_SA = np.arange(SA_binmin, SA_binmax + SA_step, SA_step)

        SA_IDs = np.empty(bins).tolist()
        for i in range(bins):
            SA_IDs[i] = []

        for i in range(bins):
            for j in range(len(SA_labels)):
                if SA_values[j] >= edges_SA[i] and SA_values[j] < edges_SA[i + 1]:
                    SAID = SA_labels[j]
                    SA_IDs[i].append(SAID)

    if 'CH4' in bin_dimensions:

        CH4_labels = np.genfromtxt(data_path + '/ch4_abs_cc_cc.txt', usecols=0,
                                   dtype=str)
        CH4_values = np.genfromtxt(data_path + '/ch4_abs_cc_cc.txt', usecols=1,
                                   dtype=float)

        CH4_binmin = 0.
        CH4_binmax = 400. # is this too low..?
        CH4_step = CH4_binmax/bins
        edges_CH4 = np.arange(CH4_binmin, CH4_binmax + CH4_step, CH4_step)

        CH4_IDs = np.empty(bins).tolist()
        for i in range(bins):
            CH4_IDs[i] = []

        for i in range(bins):
            for j in range(len(CH4_labels)):
                if CH4_values[j] >= edges_CH4[i] and CH4_values[j] < edges_CH4[i + 1]:
                    CH4ID = CH4_labels[j]
                    CH4_IDs[i].append(CH4ID)

    ###  ONE DIMENSIONAL BINNING  ###

    if len(bin_dimensions) == 1:
        freq = np.empty(bins)

        if 'CH4' in bin_dimensions:
            IDs = CH4_IDs
        if 'HV' in bin_dimensions:
            IDs = HV_IDs
        if 'SA' in bin_dimensions:
            IDs = SA_IDs

        for i in range(bins):
            freq[i] = len(IDs[i])

        mat_per_bin = np.zeros(bins)
        for i in range(bins):
            mat_per_bin[i] = freq.sum() / freq[i]

        mat_per_bin = mat_per_bin * freq.sum() / mat_per_bin.sum()

        for i in range(bins):
            mat_per_bin[i] = round(mat_per_bin[i])

        weights = mat_per_bin / mat_per_bin.sum()

        w_list = weights
        ID_list = IDs

    ###  TWO DIMENSIONAL BINNING  ###

    if len(bin_dimensions) == 2:
        IDs = np.empty([bins,bins]).tolist()
        freq = np.empty([bins,bins])
        for i in range(bins):
            for j in range(bins):
                IDs[i][j] = []

        if 'CH4' in bin_dimensions and 'HV' in bin_dimensions:
            for i in range(bins):
                for j in range(bins):
                    CH4_list = CH4_IDs[i]
                    HV_list = HV_IDs[j]
                    for k in CH4_list:
                        for l in HV_list:
                            if k == l:
                                IDs[i][j].append(k)

        if 'CH4' in bin_dimensions and 'SA' in bin_dimensions:
            for i in range(bins):
                for j in range(bins):
                    CH4_list = CH4_IDs[i]
                    SA_list = SA_IDs[j]
                    for k in CH4_list:
                        for l in SA_list:
                            if k == l:
                                IDs[i][j].append(k)

        if 'SA' in bin_dimensions and 'HV' in bin_dimensions:
            for i in range(bins):
                for j in range(bins):
                    SA_list = SA_IDs[i]
                    HV_list = HV_IDs[j]
                    for k in SA_list:
                        for l in HV_list:
                            if k == l:
                                IDs[i][j].append(k)

        for i in range(bins):
            for j in range(bins):
                freq[i,j] = len(IDs[i][j])

        mat_per_bin = np.zeros([bins,bins])
        for i in range(bins):
            for j in range(bins):
                if freq[i][j] != 0.:
                    mat_per_bin[i][j] = freq.sum() / freq[i][j]

        mat_per_bin = mat_per_bin * freq.sum() / mat_per_bin.sum() 

        for i in range(bins):
            for j in range(bins):
                mat_per_bin[i][j] = round(mat_per_bin[i][j])

        weights = mat_per_bin / mat_per_bin.sum()

        w_list = []
        ID_list = []
        for i in range(bins):
            w_list = np.concatenate([w_list, weights[i,:]])
            for j in range(bins):
                ID_list = ID_list + [IDs[i][j]]

    ###  THREE DIMENSIONAL BINNING  ###

    if len(bin_dimensions) == 3:
        IDs = np.empty([bins,bins,bins]).tolist()
        freq = np.empty([bins,bins,bins])
        for i in range(bins):
            for j in range(bins):
                for k in range(bins):
                    IDs[i][j][k] = []

        if 'CH4' in bin_dimensions and 'HV' in bin_dimensions and 'SA' in bin_dimensions:
            for i in range(bins):
                for j in range(bins):
                    for k in range(bins):
                        CH4_list = CH4_IDs[i]    
                        HV_list = HV_IDs[j]
                        SA_list = SA_IDs[k]
                        for l in CH4_list:
                            for m in HV_list:
                                if l == m:
                                    for n in SA_list:
                                        if m == n:
                                            IDs[i][j][k].append(n)

        for i in range(bins):
            for j in range(bins):
                for k in range(bins):
                    freq[i,j,k] = len(IDs[i][j][k])

        weights = np.zeros([bins,bins,bins])
        for i in range(bins):
            for j in range(bins):
                for k in range(bins):
                    if freq[i][j][k] != 0.:
                        weights[i][j][k] = freq.sum() / freq[i][j][k]

        weights = weights / weights.sum() 

        w_list = []
        ID_list = []
        for i in range(bins):
            for j in range(bins):
                w_list = np.concatenate([w_list, weights[i,j,:]])
                for k in range(bins):
                    ID_list = ID_list + [IDs[i][j][k]]

    top_path = str(atom_types) + 'atmtyp_' + str(number_of_children) + 'mutants_strength' + str(int(strength*100))

    if not os.path.exists(top_path):
        os.mkdir(top_path)

#    print( len(CH4_labels) )
    print( 'making ' + str(number_of_children) + ' materials...' )

    for i in np.arange(len(CH4_labels), len(CH4_labels) + number_of_children):
        ID_bin = np.random.choice(ID_list, p=w_list)
        parent_ID = np.random.choice(ID_bin)

#        print( 'making MAT-' + str(i) + ' (parent: ' + parent_ID + ')...' )        
        
        # creating child directory...
        mat_path = top_path + '/MAT-' + str(i)
        os.mkdir(mat_path)

#        # copying parent data files...
#        shutil.copy(data_path + '/' + parent_ID + 
#                    '/force_field_mixing_rules.def', mat_path + 
#                    '/old_mixing_rules.txt')
#
        shutil.copy(data_path + '/' + parent_ID + 
                    '/pseudo_atoms.def', mat_path)  # consider charges!
        shutil.copy(data_path + '/' + parent_ID + 
                    '/force_field.def', mat_path) 
#        shutil.copy(data_path + '/' + parent_ID + '/' + parent_ID
#                    + '.cif', mat_path + '/old_cif.txt')                        
     
          

        # importing values from parent data files...
        n1, n2, ep_o, sig_o = np.genfromtxt(data_path + '/' + parent_ID + 
                                            '/force_field_mixing_rules.def',
                                            unpack=True,
                                            skip_header=7,
                                            skip_footer=9)

        cif_atyp = np.genfromtxt(data_path + '/' + parent_ID  + '/' + 
                                 parent_ID + '.cif', usecols=0,
#                                 dtype=str, skip_header=19)
                                 dtype=str, skip_header=17)

#        n1, n2, x_o, y_o, z_o, q_o = np.genfromtxt(mat_path +       ## change this for next gens
        n1, n2, x_o, y_o, z_o = np.genfromtxt(data_path + '/' + parent_ID + 
                                                   '/' + parent_ID + '.cif',
                                                   unpack=True,
#                                                   skip_header=19)  ## change this for next gens
                                                   skip_header=17)

        # opening child data files ...            
        cif_file = open(os.path.abspath(mat_path) + '/MAT-' + str(i) + 
                        '.cif', 'w')

        plog_old = data_path + '/' + parent_ID + '/parents_log.txt'
        if not os.path.isfile(plog_old):
            parents_log = open(os.path.abspath(mat_path) + '/parents_log.txt', 'w')
            parents_log.write(parent_ID + '\nMAT-' + str(i) )
        if os.path.isfile(plog_old):
            plog = np.genfromtxt(plog_old, dtype=str)
            parents_log = open(os.path.abspath(mat_path) + '/parents_log.txt', 'w')
            for f in plog:
                parents_log.write(f + '\n')
            parents_log.write('MAT-' + str(i))

        mixing_file = open(os.path.abspath(mat_path) + 
                           '/force_field_mixing_rules.def', 'w')

        # perturbing crystal lattice parameters (a,b,c)...        
        a_footer = len(x_o)+11
        b_footer = len(x_o)+10
        c_footer = len(x_o)+9

        n1, a_o = np.genfromtxt(data_path + '/' + parent_ID + '/' + parent_ID + 
                                '.cif', unpack=True,
                                skip_header=5, skip_footer=a_footer)

        n1, b_o = np.genfromtxt(data_path + '/' + parent_ID + '/' + parent_ID + 
                                '.cif', unpack=True,
                                skip_header=6, skip_footer=b_footer)

        n1, c_o = np.genfromtxt(data_path + '/' + parent_ID + '/' + parent_ID +
                                '.cif', unpack=True,
                                skip_header=7, skip_footer=c_footer)

        a_r = round(random() * (xmax - xmin) + xmin, 4)
        b_r = round(random() * (ymax - ymin) + ymin, 4)
        c_r = round(random() * (zmax - zmin) + zmin, 4)

        a_n = round(a_o + strength * (a_r - a_o), 4)
        b_n = round(b_o + strength * (b_r - b_o), 4)
        c_n = round(c_o + strength * (c_r - c_o), 4)

        # writing new crystal lattice parameters to file...                    
        cif_head = ('MAT-' + str(i) +
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
#                    '_atom_site_charge\n' +
#                    '#NET CHARGE: ???\n')
        cif_file.write(cif_head)


        # perturbing number density...
        n_o = len(x_o)
        nden_o = n_o / (a_o * b_o * c_o)

        #    nden_r = randrange(ndenmin * 10.**16, ndenmax * 10.**16) / 10.**16
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
#                print atyp
#                print atyp[]
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
                       str(atom_types + 8) + '\n' +
                       '# type interaction, parameters.    IMPORTANT:' + 
                       ' define shortest matches first, so that more ' +
                       'specific ones overwrites these\n')
        mixing_file.write(mixing_head)

        # perturbing LJ parameters (sigma, epsilon)...
        for o in range(atom_types):

            ep_r = round(random() * (epmax - epmin) + epmin, 4)
            sig_r = round(random() * (sigmax - sigmin) + sigmin, 4)

        #epmin * 10.**16, epmax * 10.**16) / 10.**16
        #        sig_r = randrange(sigmin * 10.**16, sigmax * 10.**16) / 10.**16
        #    nden_r = round(random() * (ndenmax - ndenmin) + ndenmin, 4)

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
