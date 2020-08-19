
import numpy as np
import pytest
from pytest import approx

from htsohm.mean_delta_pbc import minimum_distance_v, minimum_distance_point, max_pair_distance


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
    mdp = minimum_distance_point(np.array([0.5, 0.5, 0.05]), np.array([0.5, 0.5, -0.05]))
    assert mdp == approx((0.0, 0.0, -0.1))

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
