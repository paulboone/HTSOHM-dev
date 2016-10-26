from htsohm.material_files import random_number_density


def test_random_number_density():
    """ Verify provided number_density_limits works by passing a suitably large pair of equal
        limits and verifying that the random number is the only one possible. If the limits are
        not enforced, there will be a range, and the test will fail.
    """

    lattice_constants = {
        'a': 1,
        'b': 1,
        'c': 1,
    }
    assert random_number_density((100000,100000), lattice_constants) == 100000
