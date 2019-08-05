from random import choice, random, randrange, uniform
import pytest

from htsohm.db import Structure
from htsohm.pseudomaterial_generator.utilities import random_number_density

@pytest.fixture
def lcs():
    s = Structure()
    s.a = 1
    s.b = 1
    s.c = 1
    return s

def test_number_density_bounds(lcs):
    assert random_number_density((100000, 100000), lcs) == 100000

def test_one_site_constraint(lcs):
    assert random_number_density((0, 0), lcs) == 1
