from random import choice, random, randrange, uniform
import pytest

from htsohm.material_files import random_number_density

@pytest.fixture
def LCs():
    lattice_constants = {
        'a' : 1,
        'b' : 1,
        'c' : 1
    }
    return lattice_constants

def test_number_density_bounds():
    assert random_number_density((100000, 100000), LCs) == 100000

def test_two_site_constraint():
    assert random_number_density((0, 1), LCs) == 2 and random_number_density((0, 2), LCs) == 2
