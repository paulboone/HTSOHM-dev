
import itertools
from math import sqrt
from operator import mul

import matplotlib.pyplot as plt
import numpy as np

def lj(eps, sigma, r):
    return 4*eps * ((sigma / r) ** 12  - (sigma / r) ** 6)

lattice = 3.889
lattice = 3.

eps_a = 158.500000
sigma_a = 2.720000

eps_h = 513.264
sigma_h = 3.052

# use Lorentz-Bertholot mixing rules
eps_ha = sqrt(eps_a * eps_h)
sigma_ha = (sigma_a + sigma_h) / 2

cutoff=14

def get_radial_distances(xyz_values, cutoff):
    # gets list of radial distances to all framework atoms in one direction of 3d space
    all_atoms = itertools.product(xyz_values, repeat=3)
    all_atoms_r = [(a,np.sqrt(np.sum(np.array(a)**2))) for a in all_atoms]
    all_atoms_r.sort(key=lambda x: x[1])
    # !=0.0 is to remove references to oneself, if necessary
    all_atoms_r = list(filter(lambda x: x[1]<cutoff and x[1]!=0.0 , all_atoms_r))
    return [x[1] for x in all_atoms_r]


print("Minimum host-adsorbate energy should be near %6.4f angstroms" % (1.122 * sigma_ha))
print("Minimum adsorbate-adsorbate energy should be near %6.4f angstroms" % (1.122 * sigma_a))

print("")
## get radial distances to framework atoms
# linspace from 0.5-5.5 should be good for all lattice >= 2.5 and cutoff <= 14
ha_r = get_radial_distances(np.linspace(0.5,5.5,6) * lattice, cutoff)
ha_lj = [lj(eps_ha, sigma_ha, r) for r in ha_r]
total_ha_energy = 8*np.sum(ha_lj)
print("%d framework atoms within cutoff in one direction; %d in all directions." % (len(ha_r), 8*len(ha_r)))
print("Total framework-adsorbate energy: %f" % total_ha_energy)

## get radial distances to other adsorbates
# linspace from 1-5 should be good for all lattice >= 2.5 and cutoff <= 14
aa_r = get_radial_distances(np.linspace(0,5,6) * lattice, cutoff)
aa_lj = [lj(eps_a, sigma_a, r) for r in aa_r]
total_aa_energy = 8*np.sum(aa_lj)
print("%d adsorbate atoms within cutoff in one direction; %d in all directions." % (len(aa_r), 8*len(aa_r)))
print("Total adsorbate-adsorbate energy: %f" % (total_aa_energy))
print("Total energy: %f" % (total_ha_energy + total_aa_energy))
# /2 is because energy of pair interaction is split across pair

print("")
print("*** host-adsorbate energies")
unique_r_lj = sorted(set([(round(r,6), "%6.4f: %6.4f" % (r, ha_lj[i])) for i, r in enumerate(ha_r)]))
for _, s in unique_r_lj:
    print(s)

print("")
print("*** adsorbate-adsorbate energies")
unique_r_lj = sorted(set([(round(r,6), "%6.4f: %6.4f" % (r, aa_lj[i])) for i, r in enumerate(aa_r)]))
for _, s in unique_r_lj:
    print(s)
