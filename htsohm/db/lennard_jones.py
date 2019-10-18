from sqlalchemy import Column, ForeignKey, Integer, String, Float
from sqlalchemy.orm import relationship

from htsohm.db import Base

class LennardJones(Base):
    __tablename__ = "lennard_jones"

    id = Column(Integer, primary_key=True)
    structure_id = Column(Integer, ForeignKey("structures.id"))

    sigma = Column(Float)
    epsilon = Column(Float)

    def atom_type_index(self):
        if self.id:
            return self.id - self.structure.lennard_jones[0].id
        else:
            return self.structure.lennard_jones.index(self)

    def exclude_cols(self):
        return ['id']

    def clone(self):
        copy = super(LennardJones, self).clone()
        return copy

    def __repr__(self):
        return "(%s, sigma: %f, epsilon: %f)" % (str(self.id), self.sigma, self.epsilon)
