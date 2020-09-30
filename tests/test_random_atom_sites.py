from random import seed, random
import pytest

import numpy as np

from htsohm.generator.random import new_material
from htsohm.generator.random import random_atom_sites, find_atom_site_with_minimum_distance
from htsohm.pair_distance import minimum_distance_v, minimum_distance_point, max_pair_distance, min_pair_distance

@pytest.fixture
def config():
    return dict(
        lattice_constant_limits=[2.0,4.0],
        sigma_limits=[2.0, 6.0],
        epsilon_limits=[2.516, 342.176],
        charge_limits=[-1.0, 1.0],
        lattice_cubic=True,
        number_of_atom_types=1,
        num_atoms_limits=[4,4],
        minimum_site_distance=0.7
    )

def test_new_material__too_close_sites_fail(config):
    seed(1)
    with pytest.raises(Exception):
        m = new_material(config, attempts=1)


def test_new_material__atom_sites_are_not_too_close(config):
    seed(0)
    m = new_material(config, attempts=2)
    assert m.min_pair_distance * m.a > config['minimum_site_distance']

def test_find_atom_site_with_minimum_distance__first_site_is_ok():
    assert find_atom_site_with_minimum_distance([], distance=1.0, uc_a = 1.01) is not None

def test_find_atom_site_with_minimum_distance__too_close_sites_fail():
    assert find_atom_site_with_minimum_distance([(0.0, 0.0, 0.0)], distance=1.0, uc_a=1.01) == None

def test_find_atom_site_with_minimum_distance__found_atom_sites_are_far_enough_away():
    starter_points = [(0.0, 0.0, 0.0)]
    uc_a = 2.0
    trial_sites = [find_atom_site_with_minimum_distance(starter_points, distance=1.0, uc_a=uc_a) for _ in range(1000)]
    distances = [uc_a * min_pair_distance(starter_points + [x]) for x in trial_sites if x is not None]
    print(len(distances), distances)
    assert (np.array(distances) >= 1.0).all() == True
