# from htsohm.db.runDB_declarative import *

# standard library imports
import os
from sqlalchemy import create_engine
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker
import yaml

# Init the database
with open(os.path.join('settings', 'database.yaml'), 'r') as yaml_file:
    dbconfig = yaml.load(yaml_file)
connection_string = dbconfig['connection_string']
engine = create_engine(connection_string)
session = sessionmaker(bind=engine)()

# Import all models
from htsohm.db.base import Base
from htsohm.db.material import Material
from htsohm.db.mutation_strength import MutationStrength

# Create tables in the engine, if they don't exist already.
Base.metadata.create_all(engine)
Base.metadata.bind = engine
