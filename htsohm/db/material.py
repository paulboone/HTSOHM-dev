from itertools import product
import math
from uuid import uuid4

from sqlalchemy import ForeignKey, Column, Integer, String, Float
from sqlalchemy.orm import relationship

from htsohm.db import Base, AtomSite, AtomTypes
from htsohm.pair_distance import min_pair_distance, max_pair_distance

class Material(Base):
    __tablename__ = 'materials'

    id           = Column(Integer, primary_key=True)
    parent_id    = Column(Integer, ForeignKey('materials.id'))

    perturbation = Column(String(10))
    generation   = Column(Integer)

    a = Column(Float)
    b = Column(Float)
    c = Column(Float)


    # relationships
    parent            = relationship("Material", remote_side=[id])
    atom_sites = relationship("AtomSite", backref="material")
    atom_types = relationship("AtomTypes", backref="material")

    @classmethod
    def add_data_columns(cls, property_config):
        print(property_config)
        for pcfg in property_config:
            for field in pcfg['fields']:
                colname = pcfg['prefix'] + field
                if hasattr(cls, colname):
                    print("WARNING: attribute %s already exists on Material")
                else:
                    setattr(cls, colname, Column(colname, Float))

    def __init__(self, parent=None, a=None, b=None, c=None, atom_sites=[], atom_types=[]):
        if parent:
            self.parent = parent
            self.parent_id = parent.id

        self.a = a
        self.b = b
        self.c = c

        self.atom_sites = atom_sites
        self.atom_types = atom_types

    def clone(self):
        copy = super(Material, self).clone()

        if self.atom_types:
            for lj in self.atom_types:
                copy.atom_types.append(lj.clone())

        if self.atom_sites:
            for atom_site in self.atom_sites:
                atom_types = copy.atom_types[atom_site.atom_types.atom_type_index()]
                copy.atom_sites.append(atom_site.clone(atom_types))


        return copy

    def exclude_cols(self):
        return ['id']

    def __repr__(self):
        return "(%s: parent: %s, a: %s )" % (str(self.id), self.parent_id, self.a)

    def minimum_unit_cells(self, cutoff):
        return (math.ceil(2 * cutoff / self.a),
                math.ceil(2 * cutoff / self.b),
                math.ceil(2 * cutoff / self.c))

    @property
    def id_or_uuid(self):
        if self.id is not None:
            return self.id
        if hasattr(self, 'uuid'):
            return self.uuid
        else:
            self.uuid = uuid4()
            return self.uuid


    @property
    def volume(self):
        return self.a * self.b * self.c

    @property
    def max_pair_distance(self):
        """in fractional coordinates"""
        return max_pair_distance([a.xyz for a in self.atom_sites])

    @property
    def min_pair_distance(self):
        """in fractional coordinates"""
        return min_pair_distance([a.xyz for a in self.atom_sites])

    def min_unit_cell_a(self, min_distance):
        return min_distance / self.min_pair_distance

    @property
    def number_density(self):
        return len(self.atom_sites)/self.volume

    @property
    def total_epsilon(self):
        return sum([s.atom_types.epsilon for s in self.atom_sites])

    @property
    def epsilon_density(self):
        return self.total_epsilon / self.volume

    @property
    def allq(self):
        return [a.q for a in self.atom_sites]

    @allq.setter
    def allq(self, allq):
        for i, q in enumerate(allq):
            self.atom_sites[i].q = q
