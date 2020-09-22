from random import seed, random

import numpy as np
import pytest
from pytest import approx

from htsohm.generator.random import random_charges

def test_random_charges__total_charge_is_zero():
    seed(0)
    assert sum(random_charges(1, 1))==0.0
    assert sum(random_charges(2, 1))==0.0
    assert sum(random_charges(3, 1))==0.0
    assert sum(random_charges(10, 1))==0.0


def test_random_charges__charges_are_within_bounds():
    seed(0)
    q = np.array(random_charges(10, 1))
    assert (q <= 1).all()
    assert (q >= -1).all()

    q = np.array(random_charges(3, 1))
    assert (q <= 1).all()
    assert (q >= -1).all()
