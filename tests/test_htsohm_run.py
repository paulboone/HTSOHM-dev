import importlib
import importlib.resources
import os
import shutil

import numpy as np
import pytest
from pytest import approx

from htsohm import htsohm_run, db, load_config_file
from htsohm.db import Material, init_database
import tests

flatten = lambda l: [item for sublist in l for item in sublist]

@pytest.fixture
def use_tmp_dir(tmpdir):
    os.chdir(tmpdir)
    return tmpdir

@pytest.fixture
def gen2_vf_geo(use_tmp_dir):
    gen2_vf_geo = [0.86301378, 0.96839894, 0.85215283, 0.96551849]

    with importlib.resources.path(tests, "gen2.db") as db_path:
        shutil.copyfile(db_path, os.path.join(use_tmp_dir, "pm.db"))

    return gen2_vf_geo

@pytest.fixture
def config_path():
    with importlib.resources.path(tests, "material_config.yaml") as config_path:
        yield config_path

@pytest.fixture
def config_multid_bins_path():
    with importlib.resources.path(tests, "config_multid_bins.yaml") as config_path:
        yield config_path

@pytest.fixture
def config(config_path):
    return load_config_file(config_path)

@pytest.mark.slow
@pytest.mark.usefixtures("use_tmp_dir")
def test_htsohm_run__runs(config_path):
    ids, _, bin_materials = htsohm_run(config_path, max_generations=2, return_run_vars=True)

    # material ids–looked up via all unique bins and bin_materials mats per bin–should be the material ids
    assert set(ids) == set([1, 2, 3, 4])
    explored_bins = [i for i,v in np.ndenumerate(bin_materials) if len(v) > 0]
    mat_indices_from_bins = flatten([bin_materials[b] for b in explored_bins])
    assert set(mat_indices_from_bins) == set([0, 1, 2, 3])

@pytest.mark.slow
@pytest.mark.usefixtures("use_tmp_dir")
def test_htsohm_run__with_many_bin_dimensions_runs(config_multid_bins_path):
    ids, _, bin_materials = htsohm_run(config_multid_bins_path, max_generations=2, return_run_vars=True)

    # material ids–looked up via all unique bins and bin_materials mats per bin–should be the material ids
    assert set(ids) == set([1, 2, 3, 4])
    explored_bins = [i for i,v in np.ndenumerate(bin_materials) if len(v) > 0]
    mat_indices_from_bins = flatten([bin_materials[b] for b in explored_bins])
    assert set(mat_indices_from_bins) == set([0, 1, 2, 3])

@pytest.mark.slow
def test_htsohm_run__restarts_ok(gen2_vf_geo, config_path):
    htsohm_run(config_path, restart_generation=3, max_generations=3)

    session = db.get_session()
    assert session.query(Material).count() == 6
    vf_results = [x[0] for x in session.query(Material.vf_geo).filter(Material.id <= 4)]
    assert vf_results == approx(gen2_vf_geo, 1e-4)

@pytest.mark.slow
def test_htsohm_run__restarts_override_will_start_at_right_place(gen2_vf_geo, config_path):
    htsohm_run(config_path, restart_generation=2, max_generations=2, override_db_errors=True)

    session = db.get_session()
    assert session.query(Material).count() == 4
    vf_results = [x[0] for x in session.query(Material.vf_geo).filter(Material.id <= 4)]
    assert vf_results[0:2] == approx(gen2_vf_geo[0:2], 1e-4)
    assert vf_results[2:4] != approx(gen2_vf_geo[2:4], 1e-4)

def test_htsohm_run__restarts_override_will_delete_extra_materials(gen2_vf_geo, config_path):
    htsohm_run(config_path, restart_generation=2, override_db_errors=True, max_generations=1)

    session = db.get_session()
    vf_results = [x[0] for x in session.query(Material.vf_geo)]
    assert len(vf_results) == 2
    assert vf_results[0:2] == approx(gen2_vf_geo[0:2], 1e-4)

def test_material__with_properties_fields_creates_the_appropriate_db_columns(config):
    init_database("sqlite:///:memory:", config['properties'])
    m = Material()
    assert set(m.__table__.columns.keys()) == {'id', 'parent_id', 'perturbation', 'a', 'b', 'c',
        'generation', 'henrys_CO2', 'henrys_CO2error', 'henrys_water', 'henrys_watererror',
        'henrys_N2', 'henrys_N2error', 'vf_geo', 'vf_raspa'}
