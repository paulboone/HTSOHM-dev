
import importlib.resources
import os

import pytest

from htsohm import htsohm_run
import tests

@pytest.mark.slow
def test_htsohm_run__runs(tmpdir):
    os.chdir(tmpdir)
    with importlib.resources.path(tests, "material_config.yaml") as config_path:
        htsohm_run(config_path, max_generations=2)
