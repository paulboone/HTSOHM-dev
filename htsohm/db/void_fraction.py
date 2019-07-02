from sqlalchemy import Column, ForeignKey, Integer, String, Float

from htsohm.db import Base

class VoidFraction(Base):
    __tablename__ = "void_fractions"

    id = Column(Integer, primary_key=True)

    # relationship with `materials`
    material_id = Column(Integer, ForeignKey("materials.id"))

    # simulation input
    adsorbate      = Column(String(16))
    temperature    = Column(Float)

    # simulation output
    void_fraction = Column(Float)
    void_fraction_geo = Column(Float)
    void_fraction_zeo = Column(Float)

    def __repr__(self):
        return "(%s: %s-%s p: %s)" % (str(self.id), self.material_id,
                    self.adsorbate, self.temperature, self.void_fraction)
