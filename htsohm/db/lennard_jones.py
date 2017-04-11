
import sys
import uuid

from sqlalchemy import Column, ForeignKey, Integer, String, Float, Boolean
from sqlalchemy.sql import text
from sqlalchemy.orm import relationship

from htsohm import config
from htsohm.db import Base, session, engine

class LennardJones(Base):
    """    """
    __tablename__ = 'lennard_jones'

    id = Column(Integer, primary_key=True)
    chemical_id = Column(String(10))
    sigma = Column(Float)
    epsilon = Column(Float)

    # relationship with 'structures'
    structure_id = Column(Integer, ForeignKey('structure.id'))

    def __init__(self, chemical_id=None, sigma=None, epsilon=None):
        """   """
        self.chemical_id = chemical_id
        self.sigma = sigma
        self.epsilon = epsilon
