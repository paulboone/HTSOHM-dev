
import numpy as np
import pytest
from pytest import approx

from htsohm.mean_delta_pbc import mean_delta_pbc, minimum_distance_point, minimum_distance_v

def test_mean_delta_pbc_basics():
    points = [(0.5, 0.5, 0.4), (0.5, 0.5, 0.5), (0.5, 0.5, 0.6)]
    mean, delta = mean_delta_pbc(points)

    assert mean == approx((0.5, 0.5, 0.5))
    assert delta == approx(0.1)

def test_mean_delta_pbc_w_pbc():
    points = [(0.5, 0.5, 0.1), (0.5, 0.5, 0.0), (0.5, 0.5, 0.9)]
    mean, delta = mean_delta_pbc(points)

    assert mean == approx((0.5, 0.5, 0.0))
    assert delta == approx(0.1)

    points = [(0.5, 0.5, 0.9), (0.5, 0.5, 0.0), (0.5, 0.5, 0.1)]
    mean, delta = mean_delta_pbc(points)

    assert mean == approx((0.5, 0.5, 0.0))
    assert delta == approx(0.1)

def test_mean_delta_pbc_which_point_centered():
    points = [(0.5, 0.5, 0.5), (0.5, 0.5, 0.1), (0.5, 0.5, 0.9)]
    mean, delta = mean_delta_pbc(points)

    assert mean == approx((0.5, 0.5, 5.0/30))
    assert delta == approx(1.0/3)


def test_minimum_distance_point():
    mdp = minimum_distance_point(np.array([0.5, 0.5, 0.4]), np.array([0.5, 0.5, 0.5]))
    assert mdp == approx((0.0, 0.0, 0.1))

    mdp = minimum_distance_point(np.array([0.5, 0.5, 0.05]), np.array([0.5, 0.5, -0.05]))
    assert mdp == approx((0.0, 0.0, -0.1))

    mdp = minimum_distance_point(np.array([0.5, 0.5, 0.95]), np.array([0.5, 0.5, 0.05]))
    assert mdp == approx((0.0, 0.0, 0.1))

def test_minimum_distance_v():
    assert minimum_distance_v(0.6, 0.4) == approx(-0.2)
    assert minimum_distance_v(0.4, 0.6) == approx(0.2)

    assert minimum_distance_v(0.1, 0.9) == approx(-0.2)
    assert minimum_distance_v(0.9, 0.1) == approx(0.2)
