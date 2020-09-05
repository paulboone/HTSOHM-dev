import numpy as np
import pytest
from pytest import approx

from htsohm.db import AtomSite
from htsohm.pair_distance import min_pair_distance
from htsohm.generator.mutate import mutate_pos_to_new_pos_w_pbc, find_good_move_position

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

def test_find_good_move_position__moves_with_only_one_atom_passes():
    site = AtomSite(x=0.5,y=0.5,z=0.5)
    assert find_good_move_position(site,[], uc_a=2, mutation_strength=0.2, distance=0.7, num_trials=1) is not None


def test_find_good_move_position__moves_on_too_small_lattice_fail():
    site = AtomSite(x=0.5,y=0.5,z=0.5)
    assert find_good_move_position(site,[], uc_a=2, mutation_strength=0.2, distance=2.1, num_trials=1) is None

def test_find_good_move_position__moves_on_bcc_at_minimum_lattice_fails():
    site = AtomSite(x=0.5,y=0.5,z=0.5)
    other_sites = [AtomSite(x=0.0,y=0.0,z=0.0)]
    assert find_good_move_position(site, other_sites, mutation_strength=0.2, uc_a=1.0, distance=(3**2)/2, num_trials=1000) is None

#
def test_find_good_move_position__moves_on_bcc_distant_points_with_low_ms_maintain_distance():
    site = AtomSite(x=0.5,y=0.5,z=0.5)
    other_sites = [AtomSite(x=0.0,y=0.0,z=0.0)]
    other_positions = [x.xyz for x in other_sites]
    uc_a = 2.0
    trial_sites = [find_good_move_position(site, other_sites, mutation_strength=0.2, distance=(3**0.5)/2, uc_a=uc_a, num_trials=1) for _ in range(1000)]
    distances = [uc_a * min_pair_distance(other_positions + [x]) for x in trial_sites if x is not None]
    assert len(distances) == 1000
    assert (np.array(distances) >= (3**0.5)/2).all() == True

def test_find_good_move_position__1_trial_moves_on_bcc_distant_points_with_high_ms_maintain_distance_should_mostly_maintain_distance():
    site = AtomSite(x=0.5,y=0.5,z=0.5)
    other_sites = [AtomSite(x=0.0,y=0.0,z=0.0)]
    other_positions = [x.xyz for x in other_sites]
    uc_a = 2.0
    trial_sites = [find_good_move_position(site, other_sites, mutation_strength=1.0, distance=(3**0.5)/2, uc_a=uc_a, num_trials=1) for _ in range(1000)]
    distances = [uc_a * min_pair_distance(other_positions + [x]) for x in trial_sites if x is not None]
    assert 600 < len(distances) < 700
    assert (np.array(distances) >= (3**0.5)/2).all() == True

def test_find_good_move_position__moves_on_bcc_distant_points_with_high_ms_maintain_distance_should_maintain_distance():
    site = AtomSite(x=0.5,y=0.5,z=0.5)
    other_sites = [AtomSite(x=0.0,y=0.0,z=0.0)]
    other_positions = [x.xyz for x in other_sites]
    uc_a = 2.0
    trial_sites = [find_good_move_position(site, other_sites, mutation_strength=1.0, distance=(3**0.5)/2, uc_a=uc_a, num_trials=100) for _ in range(1000)]
    distances = [uc_a * min_pair_distance(other_positions + [x]) for x in trial_sites if x is not None]
    assert len(distances) == 1000
    assert (np.array(distances) >= (3**0.5)/2).all() == True
