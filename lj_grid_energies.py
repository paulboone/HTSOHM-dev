
import itertools
from math import sqrt
from operator import mul

import matplotlib.pyplot as plt
import numpy as np

def lj(eps, sigma, r):
    return 4*eps * ((sigma / r) ** 12  - (sigma / r) ** 6)

lattice_lims = [2.5, 7.5]
lattice = 3.45

# Host framework
eps_h = 513.264
sigmas = [1.052, 1.663, 2.884, 3.495, 4.106, 4.717, 5.327, 5.938, 6.549]

# Adsorbate eps / sigma (CH4-Methane)
eps_a = 158.500000
sigma_a = 2.720000

cutoff=14

def get_radial_distances(xyz_values, cutoff, calc_product=True):
    # gets list of radial distances to all framework atoms in one direction of 3d space
    if calc_product:
        all_atoms = itertools.product(xyz_values, repeat=3)
    else:
        all_atoms = xyz_values

    all_atoms_r = [(a,np.sqrt(np.sum(np.array(a)**2))) for a in all_atoms]
    all_atoms_r.sort(key=lambda x: x[1])
    all_atoms_r = list(filter(lambda x: x[1]<cutoff, all_atoms_r))
    return [x[1] for x in all_atoms_r]

def get_distances_energies(xyz_values, cutoff, eps, sigma, calc_product=True):
    rs = get_radial_distances(xyz_values, cutoff, calc_product)
    ljs = [lj(eps, sigma, r) for r in rs]
    total_energy = np.sum(ljs)
    return (rs, ljs, total_energy)

def print_energies(rs, ljs):
    uniques = sorted(set([(round(r,6), "%6.4f: %6.4f" % (r, ljs[i])) for i, r in enumerate(rs)]))
    for _, s in uniques:
        print(s)

def total_energy(lattice, eps_h, sigma_h, verbose=True):

    # use Lorentz-Bertholot mixing rules
    eps_ha = sqrt(eps_a * eps_h)
    sigma_ha = (sigma_a + sigma_h) / 2

    ## get radial distances to framework atoms
    # linspace from 0.5-5.5 should be good for all lattice >= 2.5 and cutoff <= 14
    ha_r, ha_lj, total_ha_energy = get_distances_energies(np.linspace(0.5,5.5,6) * lattice, cutoff, eps_ha, sigma_ha)
    total_ha_energy *= 8

    ## get radial distances to other adsorbates
    # linspace from 1-5 should be good for all lattice >= 2.5 and cutoff <= 14
    aa_r, aa_lj, aa_quadrant_energy = get_distances_energies(np.linspace(1,5,5) * lattice, cutoff, eps_a, sigma_a)
    num_adsorbate_atoms = len(aa_r)
    lateral_atoms = list(zip(np.zeros(6), np.linspace(1,6,6)*lattice))
    aa_r_lat, aa_lj_lat, aa_lat_energy = get_distances_energies(lateral_atoms, cutoff, eps_a, sigma_a, calc_product=False)

    aa_r += aa_r_lat
    aa_lj += aa_lj_lat
    num_adsorbate_atoms = len(aa_r) * 8 + len(aa_r_lat) * 6
    total_aa_energy = (aa_quadrant_energy * 8 + aa_lat_energy * 6) # / 2
    # /2 is because energy of pair interaction is split across pair

    if verbose:
        print("Minimum host-adsorbate energy should be near %6.4f angstroms" % (1.122 * sigma_ha))
        print("Minimum adsorbate-adsorbate energy should be near %6.4f angstroms" % (1.122 * sigma_a))
        print("\n%d framework atoms within cutoff in one direction; %d in all directions." % (len(ha_r), 8*len(ha_r)))
        print("Total framework-adsorbate energy: %f" % total_ha_energy)
        print("%d, %d adsorbate atoms within cutoff in one direction; %d in all directions." % (len(aa_r), len(aa_r_lat), num_adsorbate_atoms))
        print("Total adsorbate-adsorbate energy: %f" % (total_aa_energy))
        print("Total energy: %f" % (total_ha_energy + total_aa_energy))
        print("\n*** host-adsorbate energies")
        print_energies(ha_r, ha_lj)
        print("\n*** adsorbate-adsorbate energies")
        print_energies(aa_r, aa_lj)

    return total_ha_energy + total_aa_energy

# total_energy(lattice, eps_h, sigma_h)

energy_range = [-5000, 0]
lattice_a = np.linspace(lattice_lims[0], lattice_lims[1],101)

fig = plt.figure(figsize=(12,12), tight_layout=True)
ax = fig.add_subplot(1, 1, 1)
ax.set_xlim(lattice_lims[0], lattice_lims[1])
ax.set_ylim(*energy_range)
ax.set_xlabel("Lattice constant [Å]")
ax.set_ylabel("Energy [K]")
ax.set_yticks(-(energy_range[1] - energy_range[0]) * np.array(range(0,6))/5)
ax.set_yticks(-(energy_range[1] - energy_range[0]) * np.array(range(0,50 + 1))/50, minor=True)

ax.set_xticks(lattice_lims[0] + (lattice_lims[1] - lattice_lims[0]) * np.array(range(0,11))/(11 - 1))
ax.set_xticks(lattice_a, minor=True)

# if show_grid:
ax.grid(linestyle='-', color='0.9', zorder=0, which='minor')
ax.grid(linestyle='-', color='0.3', zorder=1, which='major')


ax.axhline(-295, 0, 1, lw=2, linestyle="--", color="black", label="PR energy of gas", zorder=15)

for sigma in sigmas:
    energies = [total_energy(a, eps_h, sigma, verbose=False) for a in lattice_a]
    ax.plot(lattice_a, energies, lw=2, zorder=20, label="sigma = %4.3f" % sigma)

ax.legend()

# fig.savefig("sig_eps_a_ml.png")
# plt.close(fig)
