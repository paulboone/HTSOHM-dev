
import numpy as np
import pytest
from pytest import approx

from htsohm.pair_distance import minimum_distance_v, minimum_distance_point, \
                                     max_pair_distance, min_pair_distance

def test_minimum_distance_v__within_pbcs():
    assert minimum_distance_v(0.6, 0.4) == approx(-0.2)
    assert minimum_distance_v(0.4, 0.6) == approx(0.2)

def test_minimum_distance_v__across_pbcs():
    assert minimum_distance_v(0.1, 0.9) == approx(-0.2)
    assert minimum_distance_v(0.9, 0.1) == approx(0.2)

def test_minimum_distance_point__within_pbcs():
    mdp = minimum_distance_point(np.array([0.5, 0.5, 0.4]), np.array([0.5, 0.5, 0.5]))
    assert mdp == approx((0.0, 0.0, 0.1))

def test_minimum_distance_point__across_pbcs():
    mdp = minimum_distance_point(np.array([0.5, 0.5, 0.95]), np.array([0.5, 0.5, 0.05]))
    assert mdp == approx((0.0, 0.0, 0.1))

def test_max_pair_distance__one_point_should_return_0():
    points = [(0.5, 0.5, 0.4)]
    assert max_pair_distance(points) == approx(0.0)

def test_max_pair_distance__nopbc():
    points = [(0.5, 0.5, 0.4), (0.5, 0.5, 0.5), (0.5, 0.5, 0.6)]
    assert max_pair_distance(points) == approx(0.2)

def test_max_pair_distance__w_pbc():
    points = [(0.5, 0.5, 0.1), (0.5, 0.5, 0.0), (0.5, 0.5, 0.9)]
    assert max_pair_distance(points) == approx(0.2)

def test_max_pair_distance__max_value_should_be_corner_and_center():
    """ max value is sqrt(3) / 2 """
    points = [(0.0, 0.0, 0.0), (0.5, 0.5, 0.5)]
    assert max_pair_distance(points) == approx(3**0.5 / 2)

def test_min_pair_distance__one_site_is_one():
    assert min_pair_distance([(0.5, 0.5, 0.5)]) == approx(1.0)

def test_min_pair_distance__bcc_is_sqrt3_2():
    assert min_pair_distance([(0.0, 0.0, 0.0), (0.5, 0.5, 0.5)]) == approx(3**0.5 / 2)

def test_min_pair_distance__within_pbcs():
    assert min_pair_distance([(0.4, 0.5, 0.5), (0.5, 0.5, 0.5)]) == approx(0.1)

def test_min_pair_distance__across_pbcs():
    points = [(0.1, 0.0, 0.0), (-0.1, 0.0, 0.0)]
    assert min_pair_distance(points) == approx(0.2)
