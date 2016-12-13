import os

import yaml
from sqlalchemy import Column, ForeignKey, Integer, String, Float, Boolean, PrimaryKeyConstraint

from htsohm import config
from htsohm.db import Base, session

class MutationStrength(Base):
    """Declarative class mapping to table of mutation_strengths for each bin.
    
    Attributes:
        run_id (str): identification string for run.
        generation (int): iteration in overall bin-mutate-simulate routine.
        gas_adsorption_bin (int): value representing a region of gas
            loading parameter-space.
        surface_area_bin (int): value representing a region of surface area
            parameter-space.
        void_fraction_bin (int): value representing a region of void fraction
            parameter-space.
        strength (float): value determining the degree of perturbation for
            mutating a material.

    """
    __tablename__ = 'mutation_strengths'
    # COLUMN                                                 UNITS
    run_id = Column(String(50))                            # dimm.
    generation = Column(Integer)                           # generation#
    gas_adsorption_bin = Column(Integer)                  # dimm.
    surface_area_bin = Column(Integer)                     # dimm.
    void_fraction_bin = Column(Integer)                    # dimm.
    strength = Column(Float)

    __table_args__ = (
        PrimaryKeyConstraint('run_id', 'generation', 'gas_adsorption_bin', 'surface_area_bin', 'void_fraction_bin'),
    )

    def __init__(self, run_id=None, generation=None, gas_adsorption_bin=None,
                 surface_area_bin=None, void_fraction_bin=None, strength=None):
        self.run_id = run_id
        self.generation = generation
        self.gas_adsorption_bin = gas_adsorption_bin
        self.surface_area_bin = surface_area_bin
        self.void_fraction_bin = void_fraction_bin
        self.strength = strength

    @classmethod
    def get_prior(cls, run_id, generation, gas_adsorption_bin, surface_area_bin, void_fraction_bin):
        """
        Looks for the most recent mutation_strength row. If a row doesn't exist
        for this bin, the default value is used from the configuration file.

        Args:
            cls (classmethod): here MutationStrength.__init__ .
            run_id (str): identification string for run.
            gas_adsorption_bin (int): value representing a region of gas
                loading parameter-space.
            surface_area_bin (int): value representing a region of surface area
                parameter-space.
            void_fraction_bin (int): value representing a region of void fraction
                parameter-space.

        Returns:
            ms (float): either the mutation strength specified in the mutation
                stength datatable, or the default mutation strength if there is
                no row in the datatable corresponding to the bin.

        """

        ms = session.query(MutationStrength) \
                .filter(
                    MutationStrength.run_id == run_id,
                    MutationStrength.gas_adsorption_bin == gas_adsorption_bin,
                    MutationStrength.surface_area_bin == surface_area_bin,
                    MutationStrength.void_fraction_bin == void_fraction_bin,
                    MutationStrength.generation <= generation) \
                .order_by(MutationStrength.generation.desc()) \
                .first()

        if ms:
            return ms
        else:
            return MutationStrength(run_id, generation, gas_adsorption_bin, surface_area_bin,
                                    void_fraction_bin, config['initial_mutation_strength'])
