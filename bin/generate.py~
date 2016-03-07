from random import choice, random, randrange, randint
from functools import reduce
from math import fsum
import os
import numpy as np
from math import floor

def GCD(a,b):
	while b:
		a, b = b, a % b
	return a
		
def LCM(a,b):
	return a * b // GCD(a,b)

def LCMM(*args):
	return reduce(LCM, args)
	
def generate(N, ATOM_TYPES, ndenmax=0.084086, ndenmin=0.000013907, xmax=52.392, xmin=13.098, ymax=52.392, ymin=13.098, zmax=52.392, zmin=13.098, epmax=513.264, epmin=1.258, sigmax=6.549, sigmin=1.052, qmax=0.0, elem_charge=0.0001):
#
#    if type(N) != int:
#        print( 'N must be an 
	Ntag = str(N)
	ntag = str(ndenmax)
	xtag = str(xmax)
	ytag = str(ymax)
	ztag = str(zmax)
	eptag = str(epmax)
	sigtag = str(sigmax)
	qtag = str(qmax)
    
#	top_path = (Ntag + 'mat_' + str(ATOM_TYPES) + 'atmtyp')
	top_path = 'gen0'

   
	if not os.path.exists(top_path):
		os.mkdir(top_path)
    
	#open mat_stats.txt, to track material data    
	mat_stats = open(os.path.abspath(top_path)+ '/mat_stats.txt', 'w')
	mat_stat_heading = ('\nBOUNDARIES\nNumber of particles: ' + Ntag +
                    	    '\nnumber density:   ' + ntag + '\nx-coordinate: ' +
			    xtag + '\ny-coordinate: ' + ytag + '\nz-coordinate: ' +
			    ztag + '\nEpsilon: ' + eptag + '\nSigma: ' + sigtag 
			    + '\nCharge: ' + qtag + '\n\n' +
			    '#name     number density     xdim     ydim     '+
			    'zdim     total particles     net charge\n')
	mat_stats.write(mat_stat_heading)
    
	#MAT-XXX loop...
	for i in range(N):
		mat_name = 'MAT-' + str(i)
	
	#make MAT-XXX directory
		os.mkdir(top_path+'/'+mat_name)
	
	#open .cif file
		cif_file = open(os.path.abspath(top_path) + '/'+mat_name + '/' + 
						mat_name+'.cif', 'w')
	
	#open force_field_mixing_rules.def
		mix_name = '/'+ mat_name + '/force_field_mixing_rules.def'
		mixing_rules = open(os.path.abspath(top_path) + mix_name, 'w')
	
		psu_name = '/' + mat_name + '/pseudo_atoms.def'
		pseudo_atoms = open(os.path.abspath(top_path) + psu_name, 'w')

	#open force_field.def
		ff_name = '/'+ mat_name + '/force_field.def'
		force_field = open(os.path.abspath(top_path) + ff_name, 'w')

		xdim_ = round(random() * (xmax - xmin) + xmin, 4)
		ydim_ = round(random() * (ymax - ymin) + ymin, 4)
		zdim_ = round(random() * (zmax - zmin) + zmin, 4)
		
		Nmax = int(ndenmax * xdim_ * ydim_ * zdim_)
		n_ = randrange(2, Nmax, 1)
		nden_ = round(n_ / (xdim_ * ydim_ * zdim_))
	
		cif_heading = (mat_name + 
				   '\n\nloop_\n' +
				   '_symmetry_equiv_pos_as_xyz\n' +
				   '  x,y,z\n' +
				   '_cell_length_a          ' + str(xdim_) +
				   '\n_cell_length_b          ' + str(ydim_) +
				   '\n_cell_length_c          ' + str(zdim_) + 
				   '\n_cell_angle_alpha       90.0000\n' +
				   '_cell_angle_beta        90.0000\n' +
				   '_cell_angle_gamma       90.0000\n' +
				   'loop_\n' +
				   '_atom_site_label\n' +
				   '_atom_site_type_symbol\n' +
				   '_atom_site_fract_x\n' +
				   '_atom_site_fract_y\n' +
				   '_atom_site_fract_z\n')
		cif_file.write(cif_heading)
	
		mixing_heading = ('# general rule for shifted vs truncated\n' +
                  'shifted\n' +
                  '# general rule tailcorrections\n' +
                  'no\n' +
                  '# number of defined interactions\n' +
                  str(ATOM_TYPES + 8) + '\n' +
                  '# type interaction, parameters.    IMPORTANT: define shortest matches first, so that more specific ones overwrites these\n')
		mixing_rules.write(mixing_heading)
	
		pseudo_heading = ('#number of pseudo atoms\n' + str(ATOM_TYPES + 8) + 
					  '\n#type          print    as     chem     oxidation' +
					  '     mass       charge     polarization     ' +
					  'B-factor     radii    connectivity     anisotropic' +
					  '   anisotrop-type  tinker-type\n')
		pseudo_atoms.write(pseudo_heading)

#LJ parameters
#####
#####
		ep = []
		sig = []
		q = []
	
		for i in range(ATOM_TYPES):
			epsilon = round(random() * (epmax - epmin) + epmin, 4)
			ep.append(epsilon)
			sigma = round(random() * (sigmax -sigmin) + sigmin, 4)
			sig.append(sigma)
			charge = 0
			q.append(charge)
		
		ep_ = np.asarray(ep)
		sig_ = np.asarray(sig)
		q_ = np.asarray(q)
		ID_ = np.asarray(range(0,ATOM_TYPES))
	
		ep = ep_.reshape(-1,1)
		sig = sig_.reshape(-1,1)
		q = q_.reshape(-1,1)
		ID = ID_.reshape(-1,1)

		atoms = np.hstack((ID, ep, sig, q))

		n_atoms = np.empty([0, 4])

		for i in range(n_):
			a_type = choice(range(ATOM_TYPES))
			n_atoms = np.vstack([n_atoms, atoms[a_type, :]])
	
	#a_count = np.empty([ATOM_TYPES, 1], dtype=int)
		a_count = []
		a_id = []
		for i in range(ATOM_TYPES):
			if i in n_atoms[:, 0]:
				count = list(n_atoms[:, 0]).count(i)
				if count != 0:
					a_count.append(count)
					a_id.append(i)

		if len(a_count) != 1: #more than one atom type!
    
			temp_LCM = LCMM(*a_count)

			cm_max = floor(qmax / (temp_LCM * elem_charge / min(a_count)))    # maxiumum charge multiplier

			cm_list = []                                                      
    #cm_f = -1 * sum(cm_list)

    #print atoms
			ac_len = len(a_count)

    #if ac_len == 1:
    #    print 'this will be charge 0!!!!!!!!'
    #
    #else:
			for i in range(ac_len - 1):  #randomly choose first n-1 charge multipliers
				cm_i = randint(-1 * cm_max, cm_max)
				cm_list.append(cm_i)

				atoms[a_id[i], 3] = cm_i * temp_LCM * elem_charge / a_count[i]

			cm_f = -1 * sum(cm_list)
			atoms[a_id[-1], 3] = cm_f * temp_LCM * elem_charge / a_count[-1]

			net_q = 0
			for i in range(ac_len):
				net_q = a_count[i] * atoms[a_id[i], 3] + net_q
##############################################
		else: #single atom type!
			net_q = sum(atoms[:, 3])
    
		mat_charge = str(round(net_q, 10))
#		cif_file.write('#NET CHARGE: ' + mat_charge + '\n')
		mat_X_stats = (mat_name + '     ' + str(nden_) + '     ' + str(xdim_) +
				   '     ' + str(ydim_) + '     ' + str(zdim_) + '     ' + 
				   str(n_) + '     ' + mat_charge + '\n')
		mat_stats.write(mat_X_stats)
	
		eps = n_atoms[:,1]
		sigs = n_atoms[:,2]
#		qs = n_atoms[:,3]

	#writing mixing_rules, pseudo_atoms...
		for i in range(ATOM_TYPES):
			atom_X_pseudo = ('A_' + str(int(atoms[i,0])) + '   yes   C   C   0   ' +
                     			'12.0   ' + str(atoms[i,3]) + '   0.0   0.0   ' +
                     			'1.0  1.00   0   0  absolute   0\n')
			pseudo_atoms.write(atom_X_pseudo)
		
			atom_X_mixing = ('A_' + str(int(atoms[i,0])) + ' ' +
        			             'lennard-jones ' + str(atoms[i,1]) + ' '
        			             + str(atoms[i,2]) + '\n')
			mixing_rules.write(atom_X_mixing)                    

	#writing cif...
		for i in range(n_):
			x = round(random(), 4)
			y = round(random(), 4)
			z = round(random(), 4)
		
			atom_X_cif = ('A_' + str(int(n_atoms[i,0])) + '     ' + 'C     ' + 
        			          str(x) + '     ' + str(y) + '     ' + str(z) + '\n')
			cif_file.write(atom_X_cif)

	#SUPPORTED ADSORBATES
	# name         pseudo-atoms
	# N2       :   N_n2; N_com
	# CO2      :   C_co2; O_co2
	# methane  :   CH4_sp3
	# helium   :   He
	# hydrogen :   H_h2; H_com
	# H2       :   H_h2; H_com
	
		adsorbate_mixing = ('N_n2        lennard-jones   36.0     3.31\n' +
                		    'N_com       none\n' +
       	 		     'C_co2       lennard-jones   27.0     2.80\n' +
  	                  'O_co2       lennard-jones   79.0     3.05\n' +
  	                  'CH4_sp3     lennard-jones   158.5    3.72\n' +
  	                  'He          lennard-jones   10.9     2.64\n' +
  	                  'H_h2        none\n' +
  	                  'H_com       lennard-jones   36.7     2.958\n' +
  	                  '# general mixing rule for Lennard-Jones\n' +
  	                  'Lorentz-Berthelot')
		mixing_rules.write(adsorbate_mixing)
	
		adsorbate_pseudo = ('N_n2     yes   N   N   0   14.00674   -0.4048' +
						'   0.0   1.0   0.7   0   0   relative   0\n' +
						'N_com    no    N   -   0   0.0         0.8096' +
						'   0.0   1.0   0.7   0   0   relative   0\n' +
						'C_co2    yes   C   C   0   12.0        0.70' +
						'     0.0   1.0   0.720 0   0   relative   0\n' +
						'O_co2    yes   O   O   0   15.9994    -0.35' +
						'     0.0   1.0   0.68  0   0   relative   0\n' +
						'CH4_sp3  yes   C   C   0   16.04246    0.0' +
						'      0.0   1.0   1.00  0   0   relative   0\n' +
						'He       yes   He  He  0   4.002602    0.0' +
						'      0.0   1.0   1.0   0   0   relative   0\n' +
						'H_h2     yes   H   H   0   1.00794     0.468' +
						'    0.0   1.0   0.7   0   0   relative   0\n' +
						'H_com    no    H   H   0   0.0        - 0.936' +
						'   0.0   1.0   0.7   0   0   relative   0\n')
		pseudo_atoms.write(adsorbate_pseudo)
		
		force_field_rules = ('# rules to overwrite\n0\n' +
						 '# number of defined interactions\n0\n' +
						 '# mixing rules to overwrite\n0')
		force_field.write(force_field_rules)

		cif_file.close()
		mixing_rules.close()
		pseudo_atoms.close()
		force_field.close()

