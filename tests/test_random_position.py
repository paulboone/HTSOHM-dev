import pytest
from pytest import approx

from htsohm.generator.mutate import random_position

def test_random_position__no_pbcs():
    assert random_position(0.5, 0.6, 0.1) == approx(0.51)
    assert random_position(0.5, 0.6, 0.2) == approx(0.52)
    assert random_position(0.5, 0.6, 1.0) == approx(0.60)
    assert random_position(0.5, 0.4, 0.1) == approx(0.49)
    assert random_position(0.5, 0.4, 0.2) == approx(0.48)
    assert random_position(0.5, 0.4, 1.0) == approx(0.40)

    assert random_position(0.2, 0.6, 0.1) == approx(0.24)
    assert random_position(0.2, 0.6, 0.2) == approx(0.28)
    assert random_position(0.2, 0.6, 1.0) == approx(0.60)
    assert random_position(0.2, 0.0, 0.1) == approx(0.18)
    assert random_position(0.2, 0.0, 0.2) == approx(0.16)
    assert random_position(0.2, 0.0, 1.0) == approx(0.00)

    assert random_position(0.5, 0.0, 0.1) == approx(0.45)
    assert random_position(0.5, 0.0, 0.2) == approx(0.40)
    assert random_position(0.5, 1.0, 0.1) == approx(0.55)
    assert random_position(0.5, 1.0, 0.2) == approx(0.60)

def test_random_position__should_wrap_across_boundaries():
    assert random_position(0.9, 0.1, 0.1) == approx(0.92)
    assert random_position(0.9, 0.1, 0.2) == approx(0.94)
    assert random_position(0.9, 0.1, 1.0) == approx(0.10)
    assert random_position(0.7, 0.1, 0.1) == approx(0.74)
    assert random_position(0.7, 0.1, 0.2) == approx(0.78)
    assert random_position(0.7, 0.1, 1.0) == approx(0.10)

    assert random_position(0.1, 0.9, 0.1) == approx(0.08)
    assert random_position(0.1, 0.9, 0.2) == approx(0.06)
    assert random_position(0.1, 0.9, 1.0) == approx(0.90)
    assert random_position(0.1, 0.7, 0.1) == approx(0.06)
    assert random_position(0.1, 0.7, 0.2) == approx(0.02)
    assert random_position(0.1, 0.7, 1.0) == approx(0.70)
