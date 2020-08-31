import pytest
from pytest import approx

from htsohm.generator.mutate import mutate_pos_to_new_pos_w_pbc

def test_mutate_pos_to_new_pos_w_pbc__no_pbcs():
    assert mutate_pos_to_new_pos_w_pbc(0.5, 0.6, 0.1) == approx(0.51)
    assert mutate_pos_to_new_pos_w_pbc(0.5, 0.6, 0.2) == approx(0.52)
    assert mutate_pos_to_new_pos_w_pbc(0.5, 0.6, 1.0) == approx(0.60)
    assert mutate_pos_to_new_pos_w_pbc(0.5, 0.4, 0.1) == approx(0.49)
    assert mutate_pos_to_new_pos_w_pbc(0.5, 0.4, 0.2) == approx(0.48)
    assert mutate_pos_to_new_pos_w_pbc(0.5, 0.4, 1.0) == approx(0.40)

    assert mutate_pos_to_new_pos_w_pbc(0.2, 0.6, 0.1) == approx(0.24)
    assert mutate_pos_to_new_pos_w_pbc(0.2, 0.6, 0.2) == approx(0.28)
    assert mutate_pos_to_new_pos_w_pbc(0.2, 0.6, 1.0) == approx(0.60)
    assert mutate_pos_to_new_pos_w_pbc(0.2, 0.0, 0.1) == approx(0.18)
    assert mutate_pos_to_new_pos_w_pbc(0.2, 0.0, 0.2) == approx(0.16)
    assert mutate_pos_to_new_pos_w_pbc(0.2, 0.0, 1.0) == approx(0.00)

    assert mutate_pos_to_new_pos_w_pbc(0.5, 0.0, 0.1) == approx(0.45)
    assert mutate_pos_to_new_pos_w_pbc(0.5, 0.0, 0.2) == approx(0.40)
    assert mutate_pos_to_new_pos_w_pbc(0.5, 1.0, 0.1) == approx(0.55)
    assert mutate_pos_to_new_pos_w_pbc(0.5, 1.0, 0.2) == approx(0.60)

def test_mutate_pos_to_new_pos_w_pbc__should_wrap_across_boundaries():
    assert mutate_pos_to_new_pos_w_pbc(0.9, 0.1, 0.1) == approx(0.92)
    assert mutate_pos_to_new_pos_w_pbc(0.9, 0.1, 0.2) == approx(0.94)
    assert mutate_pos_to_new_pos_w_pbc(0.9, 0.1, 1.0) == approx(0.10)
    assert mutate_pos_to_new_pos_w_pbc(0.7, 0.1, 0.1) == approx(0.74)
    assert mutate_pos_to_new_pos_w_pbc(0.7, 0.1, 0.2) == approx(0.78)
    assert mutate_pos_to_new_pos_w_pbc(0.7, 0.1, 1.0) == approx(0.10)

    assert mutate_pos_to_new_pos_w_pbc(0.1, 0.9, 0.1) == approx(0.08)
    assert mutate_pos_to_new_pos_w_pbc(0.1, 0.9, 0.2) == approx(0.06)
    assert mutate_pos_to_new_pos_w_pbc(0.1, 0.9, 1.0) == approx(0.90)
    assert mutate_pos_to_new_pos_w_pbc(0.1, 0.7, 0.1) == approx(0.06)
    assert mutate_pos_to_new_pos_w_pbc(0.1, 0.7, 0.2) == approx(0.02)
    assert mutate_pos_to_new_pos_w_pbc(0.1, 0.7, 1.0) == approx(0.70)
