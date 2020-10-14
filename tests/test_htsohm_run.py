
import importlib.resources
import os
import shutil

import pytest

from htsohm import htsohm_run
import tests

@pytest.mark.slow
def test_htsohm_run__runs(tmpdir):
    os.chdir(tmpdir)
    with importlib.resources.path(tests, "material_config.yaml") as config_path:
        htsohm_run(config_path, max_generations=2)

@pytest.mark.slow
def test_htsohm_run__restarts_ok(tmpdir):
    os.chdir(tmpdir)
    with importlib.resources.path(tests, "test_pm.db") as db_path:
        shutil.copyfile(db_path, os.path.join(tmpdir, "pm.db"))

    with importlib.resources.path(tests, "material_config.yaml") as config_path:
        htsohm_run(config_path, restart_generation=3, max_generations=3)


def test_htsohm_run__restarts_override_will_delete_extra_materials(tmpdir):
    os.chdir(tmpdir)
    with importlib.resources.path(tests, "test_pm.db") as db_path:
        shutil.copyfile(db_path, os.path.join(tmpdir, "pm.db"))

    with importlib.resources.path(tests, "material_config.yaml") as config_path:
        htsohm_run(config_path, restart_generation=2, override_db_errors=True, max_generations=1)
