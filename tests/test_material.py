import os

import pytest
import yaml

from htsohm.config import default_configuration
from htsohm.db import Material, init_database
from htsohm.db.henrys_coefficient import henrys_m_to_v

def test_henrys_m_to_v__1_mol_sites_C_in_1m3_is_0_012():
    assert henrys_m_to_v(1e30, 6.02E+23, atom_site_mass=12.0) == 0.012

def test_henrys_m_to_v__1_site_C_in_1A3_is_0_012():
    assert henrys_m_to_v(1.0, 1.0, atom_site_mass=12.0) == (1 / 6.02E+23) * 0.012  / 1e-30


def test_material__with_properties_fields_creates_the_appropriate_db_columns():
    config = default_configuration()
    config['properties'] = yaml.safe_load("""
    properties:
    - type: 'void_fraction'
      prefix: 'vf_'
      fields: ['geo']
    - type: 'henrys_coefficients'
      prefix: 'henrys_'
      adsorbates: ['CO2', 'N2', 'H2O']
      fields: ['CO2', 'CO2error', 'N2', 'N2error', 'H2O', 'H2Oerror']
    """)["properties"]

    init_database("sqlite:///pmtest.db", config['properties'])
    m = Material()
    assert set(m.__table__.columns.keys()) == {'id', 'parent_id', 'perturbation', 'a', 'b', 'c',
        'generation', 'henrys_CO2', 'henrys_CO2error', 'henrys_H2O', 'henrys_H2Oerror',
        'henrys_N2', 'henrys_N2error', 'vf_geo'}
