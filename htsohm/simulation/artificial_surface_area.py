import os

from htsohm import config
from htsohm.simulation.calculate_bin import calc_bin
from htsohm.simulation.artificial_void_fraction import pseudo_void_fraction

def pseudo_surface_area(vf, s):
     
    return -2500 * vf ** 3 * s * (vf - 1)

def run(material, simulation_config):
    vf = pseudo_void_fraction(material.structure.number_density(),
            material.structure.average_sigma())
    results = {}
    results['sa_volumetric_surface_area'] = pseudo_surface_area(
            vf, material.structure.average_sigma())
    results['surface_area_bin'] = calc_bin(
            results['sa_volumetric_surface_area'],
            *simulation_config['limits'],
            config['number_of_convergence_bins'])
    return results
