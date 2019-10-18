import math

from sqlalchemy import Column, ForeignKey, Integer, String, Float
from sqlalchemy.orm import relationship

from htsohm.db import Base
from htsohm.db.atom_sites import AtomSite
from htsohm.db.lennard_jones import LennardJones

class Structure(Base):
    __tablename__ = "structures"

    id = Column(Integer, primary_key=True)
    material_id = Column(Integer, ForeignKey("materials.id"))

    a = Column(Float)
    b = Column(Float)
    c = Column(Float)

    atom_sites = relationship("AtomSite", backref="structure")
    lennard_jones = relationship("LennardJones", backref="structure")

    def exclude_cols(self):
        return ['id']

    def __init__(self, a=None, b=None, c=None, atom_sites=[], lennard_jones=[]):
        self.a = a
        self.b = b
        self.c = c

        self.atom_sites = atom_sites
        self.lennard_jones = lennard_jones

    def clone(self):
        copy = super(Structure, self).clone()
        if self.lennard_jones:
            for lj in self.lennard_jones:
                copy.lennard_jones.append(lj.clone())

        if self.atom_sites:
            for atom_site in self.atom_sites:
                lennard_jones = copy.lennard_jones[atom_site.lennard_jones.atom_type_index()]
                copy.atom_sites.append(atom_site.clone(lennard_jones))

        return copy

    def minimum_unit_cells(self, cutoff):
        return (math.ceil(2 * cutoff / self.a),
                math.ceil(2 * cutoff / self.b),
                math.ceil(2 * cutoff / self.c))

    @property
    def volume(self):
        return self.a * self.b * self.c

    def __repr__(self):
        return "(%s: %f, %f, %f)" % (self.id, self.a, self.b, self.c)
