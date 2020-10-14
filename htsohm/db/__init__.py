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

def get_sqlite_dbcs(database_path=None):
    if database_path is None:
        dbs = glob("*.db")
        if len(dbs) == 0:
            raise FileNotFoundError("Cannot find sqlite DBCS in the current directory: %s" % os.getcwd())
        elif len(dbs) > 1:
            print("WARNING: more than one *.db file found in this directory. Using first one: %s" % dbs[0])
        database_path = dbs[0]
    return "sqlite:///%s" % database_path

def init_database(connection_string, data_config, backup=False):
    global __engine__
    global __session__

    if connection_string[0:10] == "sqlite:///":
        db_path = connection_string[10:]
        if backup and os.path.exists(db_path):
            backup_path = db_path + "." + datetime.now().isoformat() + ".backup"
            copy2(db_path, backup_path)
            print("backing up prexisting database file %s to %s" % (db_path, backup_path))

    __engine__ = create_engine(connection_string)
    __session__ = sessionmaker(bind=__engine__)()

    # Create tables in the engine, if they don't exist already.
    Material.add_data_columns(data_config)
    Base.metadata.create_all(__engine__)
    Base.metadata.bind = __engine__

    return __engine__, __session__

# Import all models
from htsohm.db.base import Base
from htsohm.db.atom_sites import AtomSite
from htsohm.db.atom_types import AtomTypes
from htsohm.db.material import Material


def delete_extra_materials(delete_after_id):
    __engine__.execute("delete from materials where id > %d" % delete_after_id)
    __engine__.execute("delete from atom_types where material_id > %d" % delete_after_id)
    __engine__.execute("delete from atom_sites where material_id > %d" % delete_after_id)
