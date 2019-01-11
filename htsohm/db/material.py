import uuid

from sqlalchemy import Column, Integer, String, Float
from sqlalchemy.orm import relationship

from htsohm.db import Base, GasLoading, SurfaceArea, VoidFraction
from htsohm.db.structure import Structure

class Material(Base):
    """Declarative class mapping to table storing material/simulation data.

    Attributes:
        id (int): database table primary_key.
        run_id (str): identification string for run.
    """
    __tablename__ = 'materials'
    
    id           = Column(Integer, primary_key=True)
    run_id       = Column(String(50))
    uuid         = Column(String(40))
    parent       = Column(String(40))

    # structure properties
    unit_cell_volume     = Column(Float)
    number_density       = Column(Float)
    average_epsilon      = Column(Float)
    average_sigma        = Column(Float)

    # relationships
    gas_loading       = relationship("GasLoading")
    surface_area      = relationship("SurfaceArea")
    void_fraction     = relationship("VoidFraction")
    structure         = relationship("Structure", uselist=False, back_populates="material")

    def __init__(self, run_id=None, parent=None, ):
        """Init material-row.

        Args:
            self (class): row in material table.
            run_id : identification string for run (default = None).

        Initializes row in materials datatable.

        """
        self.uuid = str(uuid.uuid4())
        self.parent = parent
        self.run_id = run_id
        self.structure = Structure()

    def clone(self):
        copy = super(Material, self).clone()
        return copy
