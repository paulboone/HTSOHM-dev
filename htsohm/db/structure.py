from sqlalchemy import Column, ForeignKey, Integer, String, Float
from sqlalchemy.orm import relationship

from htsohm.db import Base
from htsohm.db.atom_sites import AtomSites
from htsohm.db.lennard_jones import LennardJones

class Structure(Base):
    __tablename__ = "structure"

    id = Column(Integer, primary_key=True)
    a = Column(Float)
    b = Column(Float)
    c = Column(Float)

    # relationship with 'materials'
    material_id = Column(Integer, ForeignKey("materials.id"))
    material = relationship("Material", back_populates="structure")

    # relationship with 'atom_sites'
    atom_sites = relationship("AtomSites")

    # relationship with 'lennard_jones'
    lennard_jones = relationship("LennardJones")

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
        if self.atom_sites:
            for atom_site in self.atom_sites:
                copy.atom_sites.append(atom_site.clone())
        if self.lennard_jones:
            for lj in self.lennard_jones:
                copy.lennard_jones.append(lj.clone())
        return copy

    @property
    def volume(self):
        return self.a * self.b * self.c
        if self.lennard_jones:
            for atom_type in self.lennard_jones:
                copy.lennard_jones.append(atom_type)
        return copy


    def __repr__(self):
        return "(%d: %f, %f, %f)" % (self.id, self.a, self.b, self.c)
