import os

from sqlalchemy import create_engine
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker
import yaml

__engine__ = None
__session__ = None

def get_session():
    return __session__

def get_engine():
    return __engine__

def init_database(connection_string):
    global __engine__
    global __session__

    if 'sqlite' in connection_string:
        print(
            'WARNING: attempting to use SQLite database! Okay for local debugging\n' +
            'but will not work with multiple workers, due to lack of locking features.'
        )

    __engine__ = create_engine(connection_string)
    __session__ = sessionmaker(bind=__engine__)()

    # Create tables in the engine, if they don't exist already.
    Base.metadata.create_all(__engine__)
    Base.metadata.bind = __engine__

# Import all models
from htsohm.db.base import Base
from htsohm.db.gas_loading import GasLoading
from htsohm.db.surface_area import SurfaceArea
from htsohm.db.void_fraction import VoidFraction
from htsohm.db.material import Material
from htsohm.db.structure import Structure
from htsohm.db.atom_sites import AtomSites
from htsohm.db.lennard_jones import LennardJones
