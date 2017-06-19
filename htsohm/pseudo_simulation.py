import os

def pseudo_void_fraction(config, nd, s):
    nd_max = config['number_density_limits'][1]
    s_max = config['sigma_limits'][1]

    return 1 - nd * s / (nd_max * s_max)

def pseudo_surface_area(vf, s):
     
    return -2500 * vf ** 3 * s * (vf - 1)

def pseudo_gas_adsorption(config, vf, sa, e):
    e_max = config['epsilon_limits'][1]
    alpha = 1 - (e_max - e) / e_max
    k = 0.1

    if vf <= 0.8:
        return alpha * (375 * (vf - 0.8) + 300) + (1 - alpha) * (k * sa * (1 - vf) + vf * 35)
        
    elif vf > 0.8:
        return alpha * (-1325 * (vf - 0.8) + 300) + (1 - alpha) * (k * sa * (1 - vf) + vf * 35)
