import os

from random import choice, random, randrange, randint
from functools import reduce
from math import fsum
import numpy as np
from math import floor


def gcd(a, b):
    while b:
        a, b = b, a % b
    return a


def lcm(a, b):
    return a * b // gcd(a, b)


def lcmm(*args):
    return reduce(lcm, args)


def generate(
        N,
        ATOM_TYPES,
        RUN_ID,
        ndenmax=0.084086,
        ndenmin=0.000013907,
        xmax=52.392,
        xmin=13.098,
        ymax=52.392,
        ymin=13.098,
        zmax=52.392,
        zmin=13.098,
        epmax=513.264,
        epmin=1.258,
        sigmax=6.549,
        sigmin=1.052,
        qmax=0.0,
        elem_charge=0.0001):

    Ntag = str(N)
    ntag = str(ndenmax)
    xtag = str(xmax)
    ytag = str(ymax)
    ztag = str(zmax)
    eptag = str(epmax)
    sigtag = str(sigmax)
    qtag = str(qmax)

    wd = os.environ['HTSOHM_DIR']      # specify $HTSOHM_DIR as working directory

    run_file = open(wd + '/' + RUN_ID + '.txt',"a")
    run_file.write( "\nBOUNDARIES\n" +
                    "Number density:\t%s\t-\t%s\n" % (ndenmin, ndenmax) +
                    "Lattice constant, a:\t%s\t-\t%s\n" % (xmin, xmax) +
                    "Lattice constant, b:\t%s\t-\t%s\n" % (ymin, ymax) +
                    "Lattice constant, c:\t%s\t-\t%s\n" % (zmin, zmax) +
                    "Epsilon:\t\t%s\t-\t%s\n" % (epmin, epmax) +
                    "Sigma:\t\t\t%s\t-\t%s\n" % (sigmin, sigmax) +
                    "Elemental charge:\t%s\n" % (elem_charge))

    ff_dir = os.environ['FF_DIR']      # output force-field files to $FF_DIR
    mat_dir = os.environ['MAT_DIR']    # output .cif-files to $MAT_DIR

    for i in range(N):                           # each iteration creates a new material

        mat_name = RUN_ID + '-' + str(i)

        def_dir = ff_dir + '/' + mat_name        # directory for material's force field
        os.mkdir( def_dir )

        cif_file = open(mat_dir + "/" + mat_name + ".cif", "w")        # structure file
        mix_file = open(def_dir + '/force_field_mixing_rules.def', "w") # LJ-parameters
        for_file = open(def_dir + '/force_field.def', "w")              # to overwrite LJ
        psu_file = open(def_dir + '/pseudo_atoms.def', "w")             # define atom-types

        xdim_ = round(random() * (xmax - xmin) + xmin, 4)
        ydim_ = round(random() * (ymax - ymin) + ymin, 4)
        zdim_ = round(random() * (zmax - zmin) + zmin, 4)

        Nmax = int(ndenmax * xdim_ * ydim_ * zdim_)
        n_ = randrange(2, Nmax, 1)
        nden_ = round(n_ / (xdim_ * ydim_ * zdim_))

        cif_file.write( "\nloop_\n" +
                        "_symmetry_equiv_pos_as_xyz\n" +
                        "  x,y,z\n" +
                        "_cell_length_a\t%s\n" % (xdim_) +
                        "_cell_length_b\t%s\n" % (ydim_) +
                        "_cell_length_c\t%s\n" % (zdim_) +
                        "_cell_angle_alpha\t90.0000\n" +
                        "_cell_angle_beta\t90.0000\n" +
                        "_cell_angle_gamma\t90.0000\n" +
                        "loop_\n" +
                        "_atom_site_label\n" +
                        "_atom_site_type_symbol\n" +
                        "_atom_site_fract_x\n" +
                        "_atom_site_fract_y\n" +
                        "_atom_site_fract_z\n" )

        mix_file.write( "# general rule for shifted vs truncated\n" +
                        "shifted\n" +
                        "# general rule tailcorrections\n" +
                        "no\n" +
                        "# number of defined interactions\n" +
                        "%s\n" % (ATOM_TYPES + 8) +
                        "# type interaction, parameters.    " +
                            "IMPORTANT: define shortest matches first, so" +
                            " that more specific ones overwrites these\n" )

        psu_file.write( "#number of pseudo atoms\n" +
                        "%s\n" % (ATOM_TYPES + 8) +
                        "#type\tprint\tas\tchem\toxidation\tmass\t" +
                            "charge\tpolarization\tB-factor\tradii\t" +
                            "connectivity\tanisotropic\tanisotrop-type\t" +
                            "tinker-type\n" )

        ep = []                        # assign LJ-parameters and partial charges
        sig = []
        q = []

        for i in range(ATOM_TYPES):
            epsilon = round(random() * (epmax - epmin) + epmin, 4)
            ep.append(epsilon)
            sigma = round(random() * (sigmax - sigmin) + sigmin, 4)
            sig.append(sigma)
            charge = 0
            q.append(charge)

        ep_ = np.asarray(ep).reshape(-1, 1)
        sig_ = np.asarray(sig).reshape(-1, 1)
        q_ = np.asarray(q).reshape(-1, 1)
        ID_ = np.asarray(range(0, ATOM_TYPES)).reshape(-1, 1)

        atoms = np.hstack((ID_, ep_, sig_, q_))

        n_atoms = np.empty([0, 4])

        for i in range(n_):
            a_type = choice(range(ATOM_TYPES))
            n_atoms = np.vstack([n_atoms, atoms[a_type, :]])

        a_count = []
        a_id = []
        for i in range(ATOM_TYPES):
            if i in n_atoms[:, 0]:
                count = list(n_atoms[:, 0]).count(i)
                if count != 0:
                    a_count.append(count)
                    a_id.append(i)

        if len(a_count) != 1:  # more than one atom type!

            temp_lcm = lcmm(*a_count)

            # maxiumum charge multiplier
            cm_max = floor(qmax / (temp_lcm * elem_charge / min(a_count)))

            cm_list = []
    #cm_f = -1 * sum(cm_list)

            ac_len = len(a_count)

    # if ac_len == 1:
    #    print 'this will be charge 0!!!!!!!!'
    #
    # else:
            for i in range(
                    ac_len -
                    1):  # randomly choose first n-1 charge multipliers
                cm_i = randint(-1 * cm_max, cm_max)
                cm_list.append(cm_i)

                atoms[a_id[i], 3] = cm_i * temp_lcm * elem_charge / a_count[i]

            cm_f = -1 * sum(cm_list)
            atoms[a_id[-1], 3] = cm_f * temp_lcm * elem_charge / a_count[-1]

            net_q = 0
            for i in range(ac_len):
                net_q = a_count[i] * atoms[a_id[i], 3] + net_q

        else:  # single atom type!
            net_q = sum(atoms[:, 3])

        mat_charge = str(round(net_q, 10))

        eps = n_atoms[:, 1]
        sigs = n_atoms[:, 2]
#        qs = n_atoms[:,3]

    # writing mixing_rules, pseudo_atoms...
        for i in range(ATOM_TYPES):
            psu_file.write( "A_%s\tyes\tC\tC\t0\t12.\t" % (int(atoms[i,0])) +
                                "%s\t0\t0\t1\t1\t0\t0\t" % (atoms[i,3]) +
                                "absolute\t0\n" )

            mix_file.write( "A_%s\tlennard-jones\t" % (int(atoms[i,0])) +
                                "%s\t%s\n" % (atoms[i,1], atoms[i,2]) )

    # writing cif...
        for i in range(n_):
            x = round(random(), 4)
            y = round(random(), 4)
            z = round(random(), 4)

            cif_file.write( "A_%s\tC\t" % (int(n_atoms[i,0])) +
                                "%s\t%s\t%s\n" % (x, y, z) )

# SUPPORTED ADSORBATES
# name         pseudo-atoms
# N2       :   N_n2; N_com
# CO2      :   C_co2; O_co2
# methane  :   CH4_sp3
# helium   :   He
# hydrogen :   H_h2; H_com
# H2       :   H_h2; H_com

        mix_file.write( "N_n2\tlennard-jones\t36.0\t3.31\n" +
                        "N_com\tnone\n" +
                        "C_co2\tlennard-jones\t27.0\t2.80\n" +
                        "O_co2\tlennard-jones\t79.0\t3.05\n" +
                        "CH4_sp3\tlennard-jones\t158.5\t3.72\n" +
                        "He\tlennard-jones\t10.9\t2.64\n" +
                        "H_h2\tnone\n" +
                        "H_com\tlennard-jones\t36.7\t2.958\n" +
                        "# general mixing rule for Lennard-Jones\n" +
                        "Lorentz-Berthelot" )

        psu_file.write( "N_n2\tyes\tN\tN\t0\t14.00674\t-0.4048\t0.0\t1.0\t" +
                            "0.7\t0\t0\trelative\t0\n" +
                        "N_com\tno\tN\t-\t0\t0.0\t0.8096\t0.0\t1.0\t0.7\t" +
                            "0\t0\trelative\t0\n" +
                        "C_co2\tyes\tC\tC\t0\t12.0\t0.70\t0.0\t1.0\t0.720" +
                            "\t0\t0\trelative\t0\n" +
                        "O_co2\tyes\tO\tO\t0\t15.9994\t-0.35\t0.0\t1.0\t" +
                            "0.68\t0\t0\trelative\t0\n" +
                        "CH4_sp3\tyes\tC\tC\t0\t16.04246\t0.0\t0.0\t1.0\t" +
                            "1.00\t0\t0\trelative\t0\n" +
                        "He\tyes\tHe\tHe\t0\t4.002602\t0.0\t0.0\t1.0\t1.0\t" +
                            "0\t0\trelative\t0\n" +
                        "H_h2\tyes\tH\tH\t0\t1.00794\t0.468\t0.0\t1.0\t0.7" +
                            "\t0\t0\trelative\t0\n" +
                        "H_com\tno\tH\tH\t0\t0.0\t-0.936\t0.0\t1.0\t0.7\t0" +
                            "\t0\trelative\t0\n" )

        for_file.write( "# rules to overwrite\n" +
                        "0\n" +
                        "# number of defined interactions\n" +
                        "0\n" +
                        "# mixing rules to overwrite\n" +
                        "0" )
#        force_field.write(force_field_rules)

        cif_file.close()
        mix_file.close()
        psu_file.close()
        for_file.close()
