from sqlalchemy import Column, ForeignKey, Integer, String, Float
from sqlalchemy.orm import relationship

from htsohm.db import Base

class AtomSite(Base):
    __tablename__ = "atom_sites"

    id = Column(Integer, primary_key=True)
    material_id = Column(Integer, ForeignKey("materials.id"), index=True)
    atom_types_id = Column(Integer, ForeignKey("atom_types.id"))

    x = Column(Float)
    y = Column(Float)
    z = Column(Float)
    q = Column(Float)

    atom_types = relationship("AtomTypes")

    def exclude_cols(self):
        return ['id']

    def clone(self, atom_types):
        copy = super(AtomSite, self).clone()
        copy.atom_types = atom_types
        return copy

    def __repr__(self):
        return "(%s, %f, %f, %f, %f)" % (str(self.id), self.x, self.y, self.z, self.q)

    @property
    def xyz(self):
        return (self.x, self.y, self.z)

    @xyz.setter
    def xyz(self, xyz):
        self.x, self.y, self.z = xyz
