
import sys
import uuid

from sqlalchemy import Column, ForeignKey, Integer, String, Float, Boolean
from sqlalchemy.sql import text

from htsohm.db import Base, session, engine

class Material(Base):
    __tablename__ = 'materials'
    # COLUMN                                                 UNITS
    id = Column(Integer, primary_key=True)                 # dimm.
    run_id = Column(String(50))                            # dimm.
    uuid = Column(String(40))
    parent_id = Column(Integer)                            # dimm.
    generation = Column(Integer)                           # generation#
    generation_index = Column(Integer)                     # index order of row in generation

    # retest columns
    retest_num = Column(Integer, default=0)
    retest_methane_loading_sum = Column(Float, default=0)
    retest_surface_area_sum = Column(Float, default=0)
    retest_void_fraction_sum = Column(Float, default=0)
    retest_passed = Column(Boolean)                        # will be NULL if retest hasn't been run

    # data collected
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

    # bins
    methane_loading_bin = Column(Integer)                  # dimm.
    surface_area_bin = Column(Integer)                     # dimm.
    void_fraction_bin = Column(Integer)                    # dimm.


    def __init__(self, run_id=None, ):
        self.uuid = str(uuid.uuid4())
        self.run_id = run_id

    @property
    def bin(self):
        return [self.methane_loading_bin, self.surface_area_bin, self.void_fraction_bin]

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

        rows = engine.connect().execute(
            sql,
            gen=self.generation,
            ml_bin=self.methane_loading_bin,
            sa_bin=self.surface_area_bin,
            vf_bin=self.void_fraction_bin
        ).fetchall()

        return len([ r for r in rows if r.in_bin ]) / len(rows)

    def calculate_retest_result(self, tolerance):
        ml_o = self.absolute_volumetric_loading    # initally-calculated values
        sa_o = self.volumetric_surface_area
        vf_o = self.helium_void_fraction

        ml_mean = self.retest_methane_loading_sum / self.retest_num
        sa_mean = self.retest_surface_area_sum / self.retest_num
        vf_mean = self.retest_void_fraction_sum / self.retest_num

        retest_failed = (abs(ml_mean - ml_o) >= tolerance * ml_o or
                         abs(sa_mean - sa_o) >= tolerance * sa_o or
                         abs(vf_mean - vf_o) >= tolerance * vf_o)

        return not retest_failed
