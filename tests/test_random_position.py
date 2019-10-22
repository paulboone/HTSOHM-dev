import pytest
from pytest import approx

from htsohm.generator.mutate import random_position

def test_random_position():
    assert random_position(0.5, 0.6, 0.1) == 0.51
    assert random_position(0.5, 0.6, 0.2) == 0.52

def test_random_position_boundaries():
    assert random_position(0.5, 1.0, 0.1) == 0.45
    assert random_position(0.5, 1.0, 0.2) == 0.40
    assert random_position(0.9, 0.1, 0.1) == 0.92
    assert random_position(0.9, 0.3, 0.5) == approx(0.10, 0.01)
