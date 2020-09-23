from random import seed, randrange

import numpy as np
import pytest
from pytest import approx

from htsohm.generator.random import random_charges
from htsohm.generator.mutate import mutate_charge, mutate_charges

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


def test_mutate_charge__total_charge_is_zero():
    seed(0)
    ms=0.2
    qlimits = (-1.0, 1.0)
    assert sum(mutate_charge([0.0], 0, ms, qlimits)) == 0.0
    assert sum(mutate_charge([-1.0, 1.0], 0, ms, qlimits)) == 0.0

    seed(0)
    assert sum(mutate_charge([-0.5, -0.5, 1.0], 0, ms, qlimits)) == 0.0
    assert sum(mutate_charge([0.5, 0.5, -1.0], 0, ms, qlimits)) == 0.0

def test_mutate_charge__charge_index_charge_is_changed():
    seed(0)
    ms=0.2
    qlimits = (-1.0, 1.0)
    assert mutate_charge([0.0, 0.0], 0, ms, qlimits)[0] != 0.0
    assert mutate_charge([0.5, -0.5], 0, ms, qlimits)[0] != 0.5

def test_mutate_charge__charge_index_charge_is_within_bounds():
    seed(0)
    ms=1.0
    qlimits = (-1.0, 1.0)
    assert mutate_charge([1.0, -1.0], 0, ms, qlimits)[0] <= 1.0
    assert mutate_charge([1.0, -1.0], 1, ms, qlimits)[0] >= -1.0
    assert mutate_charge([0.5, -0.5], 0, ms, qlimits)[0] <= 1.0

def test_random_and_mutate_charges_aggregate():
    ms=1.0
    qlimits = (-1.0, 1.0)
    for _ in range(1000):
        sd = randrange(10000000)
        seed(sd)
        q = np.array(random_charges(randrange(6) + 2, qlimits[1]))
        assert sum(q) == approx(0.0)
        assert (q <= 1).all()
        assert (q >= -1).all()
        qn = np.array(mutate_charges(q, ms, qlimits))
        assert sum(qn) == approx(0.0)
        assert (qn <= 1).all()
        assert (qn >= -1).all()
