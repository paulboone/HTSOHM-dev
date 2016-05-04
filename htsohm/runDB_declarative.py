import os
import sys
import yaml

from sqlalchemy import Column, ForeignKey, Integer, String, Float
from sqlalchemy import create_engine
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker

Base = declarative_base()

class RunData(Base):
    __tablename__ = 'data_table'
    # COLUMN                                                 UNITS
    id = Column(Integer, primary_key=True)                 # dimm.
    run_id = Column(String(50))                            # dimm.
    material_id = Column(Integer)                          # dimm.
    generation = Column(Integer)                           # generation#
    absolute_volumetric_loading = Column(Float)            # cm^3 / cm^3
    absolute_gravimetric_loading = Column(Float)           # cm^3 / g
    absolute_molar_loading = Column(Float)                 # mol / kg
    excess_volumetric_loading = Column(Float)              # cm^3 / cm^3
    excess_gravimetric_loading = Column(Float)             # cm^3 / g
    excess_molar_loading = Column(Float)                   # mol /kg
#    adsorbate_adsorbate_desorption = Column(Float)         # kJ / mol
#    host_adsorbate_desorption = Column(Float)              # kJ / mol
#    total_desorption = Column(Float)                       # kJ / mol
    unit_cell_surface_area = Column(Float)                 # angstroms ^ 2
    volumetric_surface_area = Column(Float)                # m^2 / cm^3
    gravimetric_surface_area = Column(Float)               # m^2 / g
    helium_void_fraction = Column(Float)                   # dimm.
    parent_id = Column(Integer)                            # dimm.
    methane_loading_bin = Column(Integer)                  # dimm.
    surface_area_bin = Column(Integer)                     # dimm.
    void_fraction_bin = Column(Integer)                    # dimm.
    dummy_test_result = Column(String(1))                  # y = material passes
                                                           # m = parent-material fails
                                                           # n = material fails

with open('database.yaml', 'r') as yaml_file:
    dbconfig = yaml.load(yaml_file)

connection_string = dbconfig['connection_string']
engine = create_engine(connection_string)

# Create tables in the engine, if they don't exist already.
Base.metadata.create_all(engine)
Base.metadata.bind = engine

session = sessionmaker(bind=engine)()
