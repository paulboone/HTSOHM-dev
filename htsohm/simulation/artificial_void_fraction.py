import os

from htsohm import config
from htsohm.simulation.calculate_bin import calc_bin

def pseudo_void_fraction(nd, s):
    nd_max = config['number_density_limits'][1]
    s_max = config['sigma_limits'][1]

    return 1 - nd * s / (nd_max * s_max)

def run(material, simulation_config):
    print(simulation_config)
    results = {}
    results['vf_helium_void_fraction'] = pseudo_void_fraction(
            material.structure.number_density(),material.structure.average_sigma())
    results['void_fraction_bin'] = calc_bin(
            results['vf_helium_void_fraction'],
            *simulation_config['limits'],
            config['number_of_convergence_bins'])
    return results
