from sqlalchemy import Column, ForeignKey, Integer, String, Float
from sqlalchemy.orm import relationship

from htsohm.db import Base

class AtomSites(Base):
    __tablename__ = "atom_sites"

    id = Column(Integer, primary_key=True)
    atom_type = Column(String(10))
    x = Column(Float)
    y = Column(Float)
    z = Column(Float)
    q = Column(Float)

    structure_id = Column(Integer, ForeignKey("structure.id"))
