import os

from htsohm import config
from htsohm.simulation.calculate_bin import calc_bin
from htsohm.simulation.artificial_void_fraction import pseudo_void_fraction
from htsohm.simulation.artificial_surface_area import pseudo_surface_area

def pseudo_gas_adsorption(vf, sa, e):
    e_max = config['epsilon_limits'][1]
    alpha = 1 - (e_max - e) / e_max
    k = 0.1

    if vf <= 0.8:
        return alpha * (375 * (vf - 0.8) + 300) + (1 - alpha) * (k * sa * (1 - vf) + vf * 35)
        
    elif vf > 0.8:
        return alpha * (-1325 * (vf - 0.8) + 300) + (1 - alpha) * (k * sa * (1 - vf) + vf * 35)

def run(material, simulation_config):
    vf = pseudo_void_fraction(material.structure.number_density(),
            material.structure.average_sigma())
    sa = pseudo_surface_area(vf, material.structure.average_sigma())

    results = {}
    results['ga0_absolute_volumetric_loading'] = pseudo_gas_adsorption(
            vf, sa, material.structure.average_epsilon())
    results['gas_adsorption_bin'] = calc_bin(
            results['ga0_absolute_volumetric_loading'],
            *simulation_config['limits'],
            config['number_of_convergence_bins'])
    return results
