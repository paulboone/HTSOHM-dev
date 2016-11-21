from random import choice, random, randrange, uniform
import pytest

from htsohm.material_files import random_number_density

@pytest.fixture
def LCs():
    #==============================================================================
    # default limits (see htsohm.sample.yaml)
    LC_limits = [25.6, 51.2]

    #==============================================================================
    # pick random lattice-constans within limits
    LCs = {
        'a'           : uniform(*LC_limits),
        'b'           : uniform(*LC_limits),
        'c'           : uniform(*LC_limits)
    }
    LCs['volume'] = LCs['a'] * LCs['b'] * LCs['c']
    return LCs

@pytest.fixture
def ND_limits():
    #==========================================================================
    # default limits (see htsohm.sample.yaml)
    limits = [1.49e-05, 0.02122]
    return limits

def test_boundaries(LCs, ND_limits):
#    for i in range(0, 10):
     number_density = random_number_density(ND_limits, LCs) / LCs['volume']
     assert number_density >= ND_limits[0] and number_density <= ND_limits[1]

def test_two_site_constraint(LCs, ND_limits):
#    for i in range(0, 10):
    assert random_number_density(ND_limits, LCs) >= 2
