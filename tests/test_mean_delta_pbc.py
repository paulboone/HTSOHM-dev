
from itertools import permutations

import numpy as np
import pytest
from pytest import approx

from htsohm.mean_delta_pbc import center_delta_pbc, minimum_distance_point, minimum_distance_v

def test_center_delta_pbc__nopbc():
    points = [(0.5, 0.5, 0.4), (0.5, 0.5, 0.5), (0.5, 0.5, 0.6)]
    center, delta = center_delta_pbc(points)

    assert center == approx((0.5, 0.5, 0.5))
    assert delta == approx(0.1)

def test_center_delta_pbc__w_pbc_point_order_should_not_matter():
    original_points = [(0.5, 0.5, 0.1), (0.5, 0.5, 0.0), (0.5, 0.5, 0.9)]

    for points in permutations(original_points):
        center, delta = center_delta_pbc(points)

        assert center == approx((0.5, 0.5, 0.0))
        assert delta == approx(0.1)

def test_center_delta_pbc__center_is_smallest_sphere_of_radius_delta():

    points = [(0.0, 0.0, 0.0), (0.0, 0.0, 0.01), (0.5, 0.5, 0.5)]
    center, delta = center_delta_pbc(points)

    assert center == approx((0.25, 0.25, 0.25))
    assert delta == approx(3**0.5 / 4)

    points = [(0.5, 0.5, 0.5), (0.5, 0.5, 0.1), (0.5, 0.5, 0.9)]
    center, delta = center_delta_pbc(points)

    assert center == approx((0.5, 0.5, 0.2))
    assert delta == approx(0.3)

def test_center_delta_pbc__max_value_should_be_corner_and_center():
    """ max value is sqrt(3) / 4, i.e. the largest radius sphere"""
    points = [(0.0, 0.0, 0.0), (0.5, 0.5, 0.5)]
    center, delta = center_delta_pbc(points)

    assert center == approx((0.25, 0.25, 0.25)) or center == approx((0.75, 0.75, 0.75))
    assert delta == approx(3**0.5 / 4)

def test_minimum_distance_point__within_pbcs():
    mdp = minimum_distance_point(np.array([0.5, 0.5, 0.4]), np.array([0.5, 0.5, 0.5]))
    assert mdp == approx((0.0, 0.0, 0.1))

def test_minimum_distance_point__across_pbcs():
    mdp = minimum_distance_point(np.array([0.5, 0.5, 0.05]), np.array([0.5, 0.5, -0.05]))
    assert mdp == approx((0.0, 0.0, -0.1))

    mdp = minimum_distance_point(np.array([0.5, 0.5, 0.95]), np.array([0.5, 0.5, 0.05]))
    assert mdp == approx((0.0, 0.0, 0.1))

def test_minimum_distance_v__within_pbcs():
    assert minimum_distance_v(0.6, 0.4) == approx(-0.2)
    assert minimum_distance_v(0.4, 0.6) == approx(0.2)

def test_minimum_distance_v__across_pbcs():
    assert minimum_distance_v(0.1, 0.9) == approx(-0.2)
    assert minimum_distance_v(0.9, 0.1) == approx(0.2)
