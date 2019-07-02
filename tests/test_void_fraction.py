import math

import pytest
from pytest import approx

from htsohm.void_fraction import calculate_void_fraction

def test_void_fraction_empty_crystal():
    assert calculate_void_fraction([], (4,4,4)) == 1.0

def test_void_fraction_huge_atom():
    assert calculate_void_fraction([(2,2,2,10)], (4,4,4)) == 0.0

def test_void_fraction_one_sphere():
    atoms = [(2,2,2,1)]
    assert calculate_void_fraction(atoms, (4,4,4)) == approx((4**3 - 4*math.pi/3) / 4**3, 0.01)
