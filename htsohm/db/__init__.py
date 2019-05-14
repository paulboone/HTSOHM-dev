from datetime import datetime
from glob import glob
import os
from shutil import copy2
import sys

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

def get_sqlite_dbcs():
    dbs = glob("*.db")
    if len(dbs) > 1:
        print("WARNING: more than one *.db file found in this directory. Using first one: %s" % dbs[0])
    return "sqlite:///%s" % dbs[0]

def init_database(connection_string, backup=False):
    global __engine__
    global __session__

    if connection_string[0:10] == "sqlite:///":
        print(
            'WARNING: attempting to use SQLite database! Okay for local debugging\n' +
            'but will not work with multiple workers, due to lack of locking features.',
            file=sys.stderr
        )
        db_path = connection_string[10:]
        if backup and os.path.exists(db_path):
            backup_path = db_path + "." + datetime.now().isoformat() + ".backup"
            copy2(db_path, backup_path)
            print("backing up prexisting database file %s to %s" % (db_path, backup_path))

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
