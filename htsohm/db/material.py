
import sys
import uuid

from sqlalchemy import Column, ForeignKey, Integer, String, Float, Boolean
from sqlalchemy.sql import text

from htsohm import config
from htsohm.db import Base, session, engine

class Material(Base):
    """Declarative class mapping to table storing material/simulation data.

    Attributes:
        id (int): database table primary_key.
        run_id (str): identification string for run.
        uuid (str): unique identification string for material.
        parent_id (int): uuid of parent mutated to create material.
        generation (int): iteration in overall bin-mutate-simulate routine.
        generation_index (int): order material was created in generation (used
            to determine when all materials appear in database for a particular
            generation).
        retest_num (int): iteration in re-test routine for statistical errors.
        retest_methane_loading_sum (float): sum of all absolute volumetric
            methane loadings calculated in re-test routine.
        retest_surface_area_sum (float): sum of all volumetric surface areas
            calculated in re-test routine.
        retest_void_fraction_sum (float): sum of all helium void fractions
            calculated in re-test routine.
        retest_passed (bool): true if the average of all re-test results is
            within the acceptable range of deviation.
        gl_absolute_volumetric_loading (float): absolute volumetric loading.
        gl_absolute_gravimetric_loading (float): absolute gravimetric loading.
        gl_absolute_molar_loading (float): absolute molar loading.
        gl_excess_volumetric_loading (float): excess volumetric loading.
        gl_excess_gravimetric_loading (float): excess gravimetric loading.
        gl_excess_molar_loading (float): excess molar loading.
        gl_host_host_avg (float): average energy of host-host interactions.
        gl_host_host_vdw (float): energy of host-host van der Waals
            interactions.
        gl_host_host_cou (float): energy of host-host electrostatic
            interactions.
        gl_adsorbate_adsorbate_avg (float): average energy of adsorbate-
            adsorbate interactions.
        gl_adsorbate_adsorbate_vdw (float): energy of adsorbate-adsorbate van
            der Waals interactions.
        gl_adsorbate_adsorbate_cou (float): energy of adsorbate-adsorbate
            electrostatic interactions.
        gl_host_adsorbate_avg (float): average energy of host-adsorbate
            interactions.
        gl_host_adsorbate_vdw (float): energy of host-adsorbate van der Waals
            interactions.
        gl_host_adsorbate_cou (float): energy of host-adsorbate electrostatic
            interactions.
        sa_unit_cell_surface_area (float): surface area of unit-cell.
        sa_volumetric_surface_area (float): surface area per unit volume.
        sa_gravimetric_surface_area (float): surface area per unit mass.
        vf_helium_void_fraction (float): void fraction measured with helium
            probe.
        methane_loading_bin (int): region of methane loading-space
            corresponding to the material's simulation results.
        surface_area_bin (int): region of surface area-space corresponding to
            the material's simulation results.
        void_fraction_bin (int): region of void fraction-space corresponding to
            the material's simulation results.

    """
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
    retest_gas_loading_sum = Column(Float, default=0)
    retest_surface_area_sum = Column(Float, default=0)
    retest_void_fraction_sum = Column(Float, default=0)
    retest_passed = Column(Boolean)                        # will be NULL if retest hasn't been run

    # data collected
    gl_absolute_volumetric_loading = Column(Float)            # cm^3 / cm^3
    gl_absolute_gravimetric_loading = Column(Float)           # cm^3 / g
    gl_absolute_molar_loading = Column(Float)                 # mol / kg
    gl_excess_volumetric_loading = Column(Float)              # cm^3 / cm^3
    gl_excess_gravimetric_loading = Column(Float)             # cm^3 / g
    gl_excess_molar_loading = Column(Float)                   # mol /kg
    gl_host_host_avg = Column(Float)                          # K
    gl_host_host_vdw = Column(Float)                          # K
    gl_host_host_cou = Column(Float)                          # K
    gl_adsorbate_adsorbate_avg = Column(Float)                # K
    gl_adsorbate_adsorbate_vdw = Column(Float)                # K
    gl_adsorbate_adsorbate_cou = Column(Float)                # K
    gl_host_adsorbate_avg = Column(Float)                     # K
    gl_host_adsorbate_vdw = Column(Float)                     # K
    gl_host_adsorbate_cou = Column(Float)                     # K
    sa_unit_cell_surface_area = Column(Float)                 # angstroms ^ 2
    sa_volumetric_surface_area = Column(Float)                # m^2 / cm^3
    sa_gravimetric_surface_area = Column(Float)               # m^2 / g
    vf_helium_void_fraction = Column(Float)                   # dimm.

    # bins
    gas_loading_bin = Column(Integer)                     # dimm.
    surface_area_bin = Column(Integer)                        # dimm.
    void_fraction_bin = Column(Integer)                       # dimm.


    def __init__(self, run_id=None, ):
        """Init material-row.

        Args:
            self (class): row in material table.
            run_id : identification string for run (default = None).

        Returns:
            Initializes row in materials datatable.

        """
        self.uuid = str(uuid.uuid4())
        self.run_id = run_id

    @property
    def bin(self):
        """Determine material's structure-property bin.

        Args:
            self (class): row in material table.

        Returns:
            The bin corresponding to a material's gas loading, void
            fraction, and surface area data and their postion in this three-
            dimension parameter-space.

        """
        return [self.gas_loading_bin, self.surface_area_bin, self.void_fraction_bin]

    def calculate_generation_index(self):
        """Determine material's generation-index.

        Args:
            self (class): row in material table.

        Returns:
            The generation-index is used to count the number of materials
            present in the database (that is to have all definition-files in
            the RASPA library and simulation data in the materials datatable).
            This attribute is used to determine when to stop adding new
            materials to one generation and start another.

        """
        return session.query(Material).filter(
                Material.run_id==self.run_id,
                Material.generation==self.generation,
                Material.id < self.id,
            ).count()

    def calculate_percent_children_in_bin(self):
        """Determine number of children in the same bin as their parent.

        Args:
            self (class): row in material table.

        Returns:
            Fraction of children in the same bin as parent (self).
        """
        sql = text("""
            select
                m.gas_loading_bin,
                m.surface_area_bin,
                m.void_fraction_bin,
                (
                    m.gas_loading_bin = p.gas_loading_bin and
                    m.surface_area_bin = p.surface_area_bin and
                    m.void_fraction_bin = p.void_fraction_bin
                ) as in_bin
            from materials m
            join materials p on (m.parent_id = p.id)
            where m.generation = :gen
              and p.gas_loading_bin = :gl_bin
              and p.surface_area_bin = :sa_bin
              and p.void_fraction_bin = :vf_bin
        """)

        rows = engine.connect().execute(
            sql,
            gen=self.generation,
            gl_bin=self.gas_loading_bin,
            sa_bin=self.surface_area_bin,
            vf_bin=self.void_fraction_bin
        ).fetchall()

        return len([ r for r in rows if r.in_bin ]) / len(rows)

    def calculate_retest_result(self, tolerance):
        """Determine if material has passed re-testing routine.

        Args:
            self (class): row in material table.
            tolerance (float): acceptable deviation as percent of originally-
                calculated value.

        Returns:
            (bool) True if material has NOT failed any of all re-tests.

        """
        simulations = config['material_properties']

        if 'gas_loading' in simulations:
            gl_o = self.gl_absolute_volumetric_loading    # initally-calculated values
            gl_mean = self.retest_gas_loading_sum / self.retest_num
        else:
            gl_o = 0
            gl_mean = 0

        if 'surface_area' in simulations:
            sa_o = self.sa_volumetric_surface_area
            sa_mean = self.retest_surface_area_sum / self.retest_num
        else:
            sa_o = 0
            sa_mean = 0

        if 'void_fraction' in simulations:
            vf_o = self.vf_helium_void_fraction
            vf_mean = self.retest_void_fraction_sum / self.retest_num
        else:
            vf_o = 0
            vf_mean = 0

        print('GL DEV\t%s' % (gl_mean - gl_o))
        print('SA DEV\t%s' % (sa_mean - sa_o))
        print('VF DEV\t%s' % (vf_mean - vf_o))

        retest_failed = (
            abs(gl_mean - gl_o) >= tolerance * gl_o and gl_o != 0 or
            abs(sa_mean - sa_o) >= tolerance * sa_o and sa_o != 0 or
            abs(vf_mean - vf_o) >= tolerance * vf_o and vf_o != 0
        )

        return not retest_failed
