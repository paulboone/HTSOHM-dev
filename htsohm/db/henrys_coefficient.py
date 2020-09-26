from sqlalchemy import Column, ForeignKey, Integer, String, Float

from htsohm.db import Base

class HenrysCoefficient(Base):
    __tablename__ = "henrys_coefficients"

    id = Column(Integer, primary_key=True)

    # relationship with `materials`
    material_id = Column(Integer, ForeignKey("materials.id"))

    # simulation output
    co2_henrys       = Column(Float)
    co2_henrys_error = Column(Float)
    h2o_henrys       = Column(Float)
    h2o_henrys_error = Column(Float)
    n2_henrys        = Column(Float)
    n2_henrys_error  = Column(Float)

    @property
    def m_to_v(self):
        st = self.material.structure
        na = 6.02E+23
        atom_site_mass = 0.012 # atomic mass of C in kg
        return (len(st.atom_sites) / na) * atom_site_mass / st.volume * 1e-30

    @property
    def co2_henrysv(self):
        return self.co2_henrys * self.m_to_v

    @property
    def co2_henrysv_error(self):
        return self.co2_henrys_error * self.m_to_v

    @property
    def h2o_henrysv(self):
        return self.h2o_henrys * self.m_to_v

    @property
    def h2o_henrysv_error(self):
        return self.h2o_henrys_error * self.m_to_v

    @property
    def n2_henrysv(self):
        return self.n2_henrys * self.m_to_v

    @property
    def n2_henrysv_error(self):
        return self.n2_henrys_error * self.m_to_v
