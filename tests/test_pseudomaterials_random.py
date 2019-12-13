from random import choice, random, randrange, uniform
import pytest

from htsohm.db import Structure

@pytest.fixture
def lcs():
    s = Structure()
    s.a = 1
    s.b = 1
    s.c = 1
    return s
