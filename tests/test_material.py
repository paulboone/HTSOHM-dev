from htsohm.simulation.simulate.henrys_coefficients import henrys_m_to_v

def test_henrys_m_to_v__1_mol_sites_C_in_1m3_is_0_012():
    assert henrys_m_to_v(1e30, 6.02E+23, atom_site_mass=12.0) == 0.012

def test_henrys_m_to_v__1_site_C_in_1A3_is_0_012():
    assert henrys_m_to_v(1.0, 1.0, atom_site_mass=12.0) == (1 / 6.02E+23) * 0.012  / 1e-30
