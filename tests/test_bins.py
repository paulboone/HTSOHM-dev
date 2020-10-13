from htsohm.bins import calc_bin, calc_bins

def test_calc_bin__values_within_bounds_are_binned_correctly():
    assert calc_bin(0.0, bound_min=0, bound_max=10, bins=10) == 0
    assert calc_bin(0.5, bound_min=0, bound_max=10, bins=10) == 0
    assert calc_bin(3.00, bound_min=0, bound_max=10, bins=10) == 3
    assert calc_bin(4.99, bound_min=0, bound_max=10, bins=10) == 4
    assert calc_bin(9.5, bound_min=0, bound_max=10, bins=10) == 9
    assert calc_bin(10.0, bound_min=0, bound_max=10, bins=10) == 9

def test_calc_bin__value_below_min_is_0():
    assert calc_bin(-1, bound_min=0, bound_max=10, bins=10) == 0

def test_calc_bin__value_above_max_is_last_bin():
    assert calc_bin(4000, bound_min=0, bound_max=10, bins=10) == 9

def test_calc_bins__tuples_are_binned_correctly():
    ranges = [(0,10),(0,1),(0,800)]
    assert calc_bins([(0,0,0), (4, 0.7001, 650), (10, 1, 800)], 10, ranges) == [(0,0,0), (4,7,8), (9,9,9)]
