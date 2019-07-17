import math

from sqlalchemy import Column, ForeignKey, Integer, String, Float
from sqlalchemy.orm import relationship

from htsohm.db import Base
from htsohm.db.atom_sites import AtomSite
from htsohm.db.lennard_jones import LennardJones

class Structure(Base):
    __tablename__ = "structures"

    id = Column(Integer, primary_key=True)
    a = Column(Float)
    b = Column(Float)
    c = Column(Float)

    # relationship with 'materials'
    material_id = Column(Integer, ForeignKey("materials.id"))
    material = relationship("Material", back_populates="structure")

    # relationship with 'atom_sites'
    atom_sites = relationship("AtomSite")

    # relationship with 'lennard_jones'
    lennard_jones = relationship("LennardJones")

    def get_lennard_jones(self, atom_type):
        a = [lj for lj in self.lennard_jones if lj.atom_type == atom_type][0]
        return a

    def map_atom_sites_to_lj(self):
        # assign lennard jones by id so relationship will work
        for a in self.atom_sites:
            a.lennard_jones = self.get_lennard_jones(a.atom_type)

    def exclude_cols(self):
        return ['id']

    def __init__(self, a=None, b=None, c=None, atom_sites=[], lennard_jones=[]):
        self.a = a
        self.b = b
        self.c = c

        self.atom_sites = atom_sites
        self.lennard_jones = lennard_jones

        self.map_atom_sites_to_lj()

    def clone(self):
        copy = super(Structure, self).clone()
        if self.lennard_jones:
            for lj in self.lennard_jones:
                copy.lennard_jones.append(lj.clone())

        if self.atom_sites:
            for atom_site in self.atom_sites:
                copy.atom_sites.append(atom_site.clone())

        copy.map_atom_sites_to_lj()
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
