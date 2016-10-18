import os

import yaml
from sqlalchemy import Column, ForeignKey, Integer, String, Float, Boolean, PrimaryKeyConstraint

from htsohm import config
from htsohm.db import Base, session

class MutationStrength(Base):
    __tablename__ = 'mutation_strengths'
    # COLUMN                                                 UNITS
    run_id = Column(String(50))                            # dimm.
    generation = Column(Integer)                           # generation#
    methane_loading_bin = Column(Integer)                  # dimm.
    surface_area_bin = Column(Integer)                     # dimm.
    void_fraction_bin = Column(Integer)                    # dimm.
    strength = Column(Float)

    __table_args__ = (
        PrimaryKeyConstraint('run_id', 'generation', 'methane_loading_bin', 'surface_area_bin', 'void_fraction_bin'),
    )

    def __init__(self, run_id=None, generation=None, methane_loading_bin=None,
                 surface_area_bin=None, void_fraction_bin=None, strength=None):
        self.run_id = run_id
        self.generation = generation
        self.methane_loading_bin = methane_loading_bin
        self.surface_area_bin = surface_area_bin
        self.void_fraction_bin = void_fraction_bin
        self.strength = strength

    @classmethod
    def get_prior(cls, run_id, generation, methane_loading_bin, surface_area_bin, void_fraction_bin):
        """
        Looks for the most recent mutation_strength row for the passed run_id, generation, and bins
        and returns it. If a row doesn't exist for this bin, return a
        """

        ms = session.query(MutationStrength) \
                .filter(
                    MutationStrength.run_id == run_id,
                    MutationStrength.methane_loading_bin == methane_loading_bin,
                    MutationStrength.surface_area_bin == surface_area_bin,
                    MutationStrength.void_fraction_bin == void_fraction_bin,
                    MutationStrength.generation <= generation) \
                .order_by(MutationStrength.generation.desc()) \
                .first()

        if ms:
            return ms
        else:
            return MutationStrength(run_id, generation, methane_loading_bin, surface_area_bin,
                                    void_fraction_bin, config['initial_mutation_strength'])
