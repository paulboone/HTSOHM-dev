from sqlalchemy import Column, ForeignKey, Integer, String, Float
from sqlalchemy.orm import relationship

from htsohm.db import Base

class AtomTypes(Base):
    __tablename__ = "atom_types"

    id = Column(Integer, primary_key=True)
    material_id = Column(Integer, ForeignKey("materials.id"))

    sigma = Column(Float)
    epsilon = Column(Float)

    def atom_type_index(self):
        if self.id:
            return self.id - self.material.atom_types[0].id
        else:
            return self.material.atom_types.index(self)

    def exclude_cols(self):
        return ['id']

    def clone(self):
        copy = super(AtomTypes, self).clone()
        return copy

    def __repr__(self):
        return "(%s, sigma: %f, epsilon: %f)" % (str(self.id), self.sigma, self.epsilon)
