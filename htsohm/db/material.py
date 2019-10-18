import uuid

from sqlalchemy import ForeignKey, Column, Integer, String, Float
from sqlalchemy.orm import relationship

from htsohm.db import Base, GasLoading, SurfaceArea, VoidFraction, AtomSite, LennardJones
from htsohm.db.structure import Structure

class Material(Base):
    __tablename__ = 'materials'

    id           = Column(Integer, primary_key=True)
    uuid         = Column(String(40))
    parent_id    = Column(Integer, ForeignKey('materials.id'))
    perturbation = Column(String(10))
    generation   = Column(Integer)

    # structure properties
    number_density       = Column(Float)

    # relationships
    parent            = relationship("Material", remote_side=[id])
    structure         = relationship("Structure", uselist=False, backref="material", cascade="all, delete-orphan")

    gas_loading       = relationship("GasLoading", cascade="all, delete-orphan")
    surface_area      = relationship("SurfaceArea", cascade="all, delete-orphan")
    void_fraction     = relationship("VoidFraction", cascade="all, delete-orphan")

    def __init__(self, parent=None, structure=None, number_density=None):
        """Init material-row.

        Args:
            self (class): row in material table.

        Initializes row in materials datatable.

        """
        self.uuid = str(uuid.uuid4())
        if parent:
            self.parent = parent
            self.parent_id = parent.id
        if structure is None:
            self.structure = Structure()
        else:
            self.structure = structure

        if number_density:
            self.number_density = number_density

    @staticmethod
    def one_atom_new(sigma, epsilon, a, b, c):
        structure = Structure(a, b, c,
                [AtomSite(atom_type="A_0", x=1.0, y=1.0, z=1.0, q=0.0)],
                [LennardJones(atom_type="A_0", sigma=sigma, epsilon=epsilon)])
        m = Material(structure=structure)
        return m

    @staticmethod
    def cube_pore_new(sigma, epsilon, num_atoms, atom_diameter):
        # lattice constant a is calculated from number of atoms times the atom_diameter
        a = num_atoms * atom_diameter

        atom_sites = []
        for xi in range(num_atoms):
            for yi in range(num_atoms):
                for zi in range(num_atoms):
                    # only add indices that are in one of the three boundary planes, i.e index == 0
                    if min(xi, yi, zi) == 0:
                        x = xi * atom_diameter / a
                        y = yi * atom_diameter / a
                        z = zi * atom_diameter / a
                        atom_sites.append(AtomSite(atom_type="A_0", x=x, y=y, z=z, q=0.0))

        structure = Structure(a, a, a, atom_sites, [LennardJones(atom_type="A_0", sigma=sigma, epsilon=epsilon)])
        m = Material(structure=structure)
        return m



    def clone(self):
        copy = super(Material, self).clone()
        copy.parent = self
        copy.parent_id = self.id
        copy.structure = self.structure.clone()
        return copy

    def exclude_cols(self):
        return ['uuid', 'id']

    def __repr__(self):
        return "(%s: %s p: %s)" % (str(self.id), self.uuid, self.parent_id)
