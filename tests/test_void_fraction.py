import math

import pytest
from pytest import approx

from htsohm.void_fraction import calculate_void_fraction

r1volume = 4*math.pi/3

def test_void_fraction_empty_crystal():
    assert calculate_void_fraction([], (4,4,4)) == 1.0

def test_void_fraction_huge_atom():
    assert calculate_void_fraction([(2,2,2,8)], (4,4,4), points_per_angstrom=2) == 0.0

def test_void_fraction_huger_atom():
    assert calculate_void_fraction([(2,2,2,100)], (4,4,4), points_per_angstrom=2) == 0.0

def test_void_fraction_one_sphere():
    atoms = [(2,2,2,2)]
    assert calculate_void_fraction(atoms, (4,4,4)) == approx(1 - r1volume/4**3, 0.01)

def test_void_fraction_one_sphere_w_he_probe():
    atoms = [(2,2,2,2)]
    assert calculate_void_fraction(atoms, (4,4,4), probe_r=2.105) == approx(1 - r1volume*((4.105/2)/4)**3, 0.01)

def test_void_fraction_two_identical_spheres():
    atoms = [(2,2,2,2), (2,2,2,2)]
    assert calculate_void_fraction(atoms, (4,4,4)) == approx(1 - r1volume/4**3, 0.01)

def test_void_fraction_two_separate_spheres():
    atoms = [(1,1,1,2), (2.9,2.9,2.9,2)] #2.9 so it doesn't hit the right boundary
    assert calculate_void_fraction(atoms, (4,4,4)) == approx(1 - 2*r1volume/4**3, 0.01)

def test_void_fraction_one_sphere_across_pbc():
    atoms = [(0,0,0,2)]
    assert calculate_void_fraction(atoms, (4,4,4)) == approx(1 - r1volume/4**3, 0.01)

    atoms = [(4,4,4,2)]
    assert calculate_void_fraction(atoms, (4,4,4)) == approx(1 - r1volume/4**3, 0.01)
