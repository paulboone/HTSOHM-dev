from sqlalchemy import Column, ForeignKey, Integer, String, Float

from htsohm.db import Base

class GasLoading(Base):
    __tablename__ = "gas_loadings"

    id = Column(Integer, primary_key=True)

    # relationship with `materials`
    material_id = Column(Integer, ForeignKey("materials.id"))

    # simulation input
    adsorbate      = Column(String(16))
    pressure       = Column(Float)
    temperature    = Column(Float)

    # simulation output
    cycles                           = Column(Integer)
    absolute_volumetric_loading      = Column(Float)
    absolute_volumetric_loading_error = Column(Float)
