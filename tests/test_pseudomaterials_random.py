from random import choice, random, randrange, uniform
import pytest

from htsohm.db import Structure
from htsohm.pseudomaterial_generator.utilities import random_number_density

@pytest.fixture
def LCs():
    s = Structure()
    s.a = 1
    s.b = 1
    s.c = 1
    return s

def test_number_density_bounds():
    assert random_number_density((100000, 100000), LCs()) == 100000

def test_two_site_constraint():
    assert random_number_density((1, 2), LCs()) == 2

