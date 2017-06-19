
import sys
import uuid

from sqlalchemy import Column, ForeignKey, Integer, String, Float, Boolean
from sqlalchemy.sql import text
from sqlalchemy.orm import relationship

from htsohm import config
from htsohm.db import Base, session, engine
from htsohm.db.atom_sites import AtomSites
from htsohm.db.lennard_jones import LennardJones

class Structure(Base):
    """    """
    __tablename__ = 'structure'

    id = Column(Integer, primary_key=True)
    lattice_constant_a = Column(Float)
    lattice_constant_b = Column(Float)
    lattice_constant_c = Column(Float)
    
    # relationship with 'materials'
    material_id = Column(Integer, ForeignKey('materials.id'))
    material = relationship("Material", back_populates="structure")

    # relationship with 'atom_sites'
    atom_sites = relationship("AtomSites")

    # relationship with 'lennard_jones'
    lennard_jones = relationship("LennardJones")

    def __init__(self,
            lattice_constant_a=None,
            lattice_constant_b=None,
            lattice_constant_c=None,
            atom_sites=[],
            lennard_jones=[]):
        """   """
        self.lattice_constant_a = lattice_constant_a
        self.lattice_constant_b = lattice_constant_b
        self.lattice_constant_c = lattice_constant_c
        self.atom_sites = atom_sites
        self.lennard_jones = lennard_jones

    def clone(self):
        copy = super(Structure, self).clone()
        if self.atom_sites:
            for atom_site in self.atom_sites:
                copy.atom_sites.append(atom_site)
        if self.lennard_jones:
            for atom_type in self.lennard_jones:
                copy.lennard_jones.append(atom_type)
        return copy

    @property
    def volume(self):
        """Calculates unit cell volume by mutliplying lattice constants.
        """
        return self.lattice_constant_a * self.lattice_constant_b * self.lattice_constant_c

    def number_density(self):
        return len(self.atom_sites) / self.volume

    def average_sigma(self):
        sigmas = []
        for atom_type in self.lennard_jones:
            sigmas.append(atom_type.sigma)
        return sum(sigmas) / max(len(sigmas), 1)

    def average_epsilon(self):
        epsilons = []
        for atom_type in self.lennard_jones:
            epsilons.append(atom_type.epsilon)
        return sum(epsilons) / max(len(epsilons), 1)
