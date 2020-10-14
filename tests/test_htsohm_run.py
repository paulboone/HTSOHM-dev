
import importlib.resources
import os
import shutil

import pytest
from pytest import approx

from htsohm import htsohm_run, db
from htsohm.db import Material
import tests

@pytest.mark.slow
def test_htsohm_run__runs(tmpdir):
    os.chdir(tmpdir)
    with importlib.resources.path(tests, "material_config.yaml") as config_path:
        htsohm_run(config_path, max_generations=2)

@pytest.mark.slow
def test_htsohm_run__restarts_ok(tmpdir):
    os.chdir(tmpdir)

    gen2_vf_geo = [0.86301378, 0.96839894, 0.85215283, 0.96551849]

    with importlib.resources.path(tests, "gen2.db") as db_path:
        shutil.copyfile(db_path, os.path.join(tmpdir, "pm.db"))

    with importlib.resources.path(tests, "material_config.yaml") as config_path:
        htsohm_run(config_path, restart_generation=3, max_generations=3)

        session = db.get_session()
        assert session.query(Material).count() == 6
        vf_results = [x[0] for x in session.query(Material.vf_geo).filter(Material.id <= 4)]
        assert vf_results == approx(gen2_vf_geo, 1e-4)

@pytest.mark.slow
def test_htsohm_run__restarts_override_will_start_at_right_place(tmpdir):
    os.chdir(tmpdir)

    gen2_vf_geo = [0.86301378, 0.96839894, 0.85215283, 0.96551849]

    with importlib.resources.path(tests, "gen2.db") as db_path:
        shutil.copyfile(db_path, os.path.join(tmpdir, "pm.db"))

    with importlib.resources.path(tests, "material_config.yaml") as config_path:
        htsohm_run(config_path, restart_generation=2, max_generations=2, override_db_errors=True)

        session = db.get_session()
        assert session.query(Material).count() == 4
        vf_results = [x[0] for x in session.query(Material.vf_geo).filter(Material.id <= 4)]
        assert vf_results[0:2] == approx(gen2_vf_geo[0:2], 1e-4)
        assert vf_results[2:4] != approx(gen2_vf_geo[2:4], 1e-4)

def test_htsohm_run__restarts_override_will_delete_extra_materials(tmpdir):
    os.chdir(tmpdir)
    with importlib.resources.path(tests, "gen2.db") as db_path:
        shutil.copyfile(db_path, os.path.join(tmpdir, "pm.db"))

    with importlib.resources.path(tests, "material_config.yaml") as config_path:
        htsohm_run(config_path, restart_generation=2, override_db_errors=True, max_generations=1)
