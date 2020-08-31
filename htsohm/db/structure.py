import math

from sqlalchemy import Column, ForeignKey, Integer, String, Float
from sqlalchemy.orm import relationship

from htsohm.db import Base
from htsohm.db.atom_sites import AtomSite
from htsohm.db.atom_types import AtomTypes
from htsohm.max_pair_distance import max_pair_distance

class Structure(Base):
    __tablename__ = "structures"

    id = Column(Integer, primary_key=True)
    material_id = Column(Integer, ForeignKey("materials.id"))

    a = Column(Float)
    b = Column(Float)
    c = Column(Float)

    atom_sites = relationship("AtomSite", backref="structure")
    atom_types = relationship("AtomTypes", backref="structure")

    def exclude_cols(self):
        return ['id']

    def __init__(self, a=None, b=None, c=None, atom_sites=[], atom_types=[]):
        self.a = a
        self.b = b
        self.c = c

        self.atom_sites = atom_sites
        self.atom_types = atom_types

    def clone(self):
        copy = super(Structure, self).clone()
        if self.atom_types:
            for lj in self.atom_types:
                copy.atom_types.append(lj.clone())

        if self.atom_sites:
            for atom_site in self.atom_sites:
                atom_types = copy.atom_types[atom_site.atom_types.atom_type_index()]
                copy.atom_sites.append(atom_site.clone(atom_types))

        return copy

    def minimum_unit_cells(self, cutoff):
        return (math.ceil(2 * cutoff / self.a),
                math.ceil(2 * cutoff / self.b),
                math.ceil(2 * cutoff / self.c))

    @property
    def volume(self):
        return self.a * self.b * self.c

    @property
    def max_pair_distance(self):
        """in fractional coordinates"""
        return max_pair_distance([(a.x, a.y, a.z) for a in self.atom_sites])

    @property
    def number_density(self):
        return len(self.atom_sites)/self.volume

    @property
    def total_epsilon(self):
        return sum([s.atom_types.epsilon for s in self.atom_sites])

    @property
    def epsilon_density(self):
        return self.total_epsilon / self.volume

    def __repr__(self):
        return "(%s: %f, %f, %f)" % (self.id, self.a, self.b, self.c)
