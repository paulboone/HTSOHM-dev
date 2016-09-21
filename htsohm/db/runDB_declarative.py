# standard library imports
import os
import sys
import uuid

from sqlalchemy import Column, ForeignKey, Integer, String, Float, Boolean
from sqlalchemy import create_engine
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker
from sqlalchemy.sql import text
import yaml

Base = declarative_base()
class Material(Base):
    __tablename__ = 'materials'
    # COLUMN                                                 UNITS
    id = Column(Integer, primary_key=True)                 # dimm.
    run_id = Column(String(50))                            # dimm.
    uuid = Column(String(40))
    seed = Column(Boolean, default=False)
    generation = Column(Integer)                           # generation#
    generation_index = Column(Integer)                     # index order of row in generation
    absolute_volumetric_loading = Column(Float)            # cm^3 / cm^3
    absolute_gravimetric_loading = Column(Float)           # cm^3 / g
    absolute_molar_loading = Column(Float)                 # mol / kg
    excess_volumetric_loading = Column(Float)              # cm^3 / cm^3
    excess_gravimetric_loading = Column(Float)             # cm^3 / g
    excess_molar_loading = Column(Float)                   # mol /kg
    host_host_avg = Column(Float)                          # K
    host_host_vdw = Column(Float)                          # K
    host_host_cou = Column(Float)                          # K
    adsorbate_adsorbate_avg = Column(Float)                # K
    adsorbate_adsorbate_vdw = Column(Float)                # K
    adsorbate_adsorbate_cou = Column(Float)                # K
    host_adsorbate_avg = Column(Float)                     # K
    host_adsorbate_vdw = Column(Float)                     # K
    host_adsorbate_cou = Column(Float)                     # K
    unit_cell_surface_area = Column(Float)                 # angstroms ^ 2
    volumetric_surface_area = Column(Float)                # m^2 / cm^3
    gravimetric_surface_area = Column(Float)               # m^2 / g
    helium_void_fraction = Column(Float)                   # dimm.
    parent_id = Column(Integer)                            # dimm.
    methane_loading_bin = Column(Integer)                  # dimm.
    surface_area_bin = Column(Integer)                     # dimm.
    void_fraction_bin = Column(Integer)                    # dimm.
    write_check = Column(String(4))                        # 'done' flag means all .def files were written
    dummy_test_result = Column(String(4))                  # "pass" = material passes
                                                           # "fail" = material fails
    data_complete = Column(Boolean, server_default="0")    # set when all columns are populated

    def __init__(self, run_id, dummy_test_result):
        self.run_id = run_id
        self.dummy_test_result = dummy_test_result
        self.uuid = str(uuid.uuid4())

    def calculate_generation_index(self):
        return session.query(Material).filter(
                Material.run_id==self.run_id,
                Material.generation==self.generation,
                Material.id < self.id,
            ).count()

    def calculate_percent_children_in_bin(self):
        sql = text("""
            select
                m.methane_loading_bin,
                m.surface_area_bin,
                m.void_fraction_bin,
                (
                    m.methane_loading_bin = p.methane_loading_bin and
                    m.surface_area_bin = p.surface_area_bin and
                    m.void_fraction_bin = p.void_fraction_bin
                ) as in_bin
            from materials m
            join materials p on (m.parent_id = p.id)
            where m.generation = :gen
              and p.methane_loading_bin = :ml_bin
              and p.surface_area_bin = :sa_bin
              and p.void_fraction_bin = :vf_bin
        """)

        rows = session.connection.execute(
            sql,
            gen=self.generation,
            ml_bin=self.methane_loading_bin,
            sa_bin=self.surface_area_bin,
            vf_bin=self.void_fraction_bin
        ).fetchall()

        return len([ r for r in rows if r.in_bin ]) / len(rows)




with open(os.path.join('settings', 'database.yaml'), 'r') as yaml_file:
    dbconfig = yaml.load(yaml_file)
connection_string = dbconfig['connection_string']
engine = create_engine(connection_string)

# Create tables in the engine, if they don't exist already.
Base.metadata.create_all(engine)
Base.metadata.bind = engine
session = sessionmaker(bind=engine)()
