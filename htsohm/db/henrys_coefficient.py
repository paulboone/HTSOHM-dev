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
