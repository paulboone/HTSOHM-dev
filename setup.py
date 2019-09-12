#!/usr/bin/env python

from distutils.core import setup
setup(
    name = 'HTSOHM',
    version = '0.5.1',
    description = 'High-throughput Screening of Hypothetical Materials',
    author = 'Paul Boone / A Reino Kaija',
    author_email = 'paulboone@pitt.edu',
    url = 'https://github.com/paulboone/HTSOHM-dev',
    packages = ['htsohm'],
    install_requires=[
        'numpy',
        'scipy',
        'sqlalchemy',
        'click',
        'pandas',
        'pyyaml',
        'matplotlib',
    ],
    entry_points={
          'console_scripts': [
              'dps = htsohm.bin.dps:cmdline',
              'psm-prop-histogram = htsohm.bin.prop_histogram:bin_graph',
              'psm-graph-bins = htsohm.bin.bin_graph:bin_graph',
              'psm-graph-a-ml = htsohm.bin.graph_sig_eps_a_ml:bin_graph',
              'psm-output-config-files = htsohm.bin.config_files:output_config_files',
              'psm-csv = htsohm.bin.csv:output_csv',
              'psm-dof-analysis = htsohm.bin.dof_analysis:dof_analysis',
              'psm-setup-one-atom-sweep = htsohm.bin.one_atom_sweep_setup:sweep_setup',
              'psm-setup-cube-pore-sweep = htsohm.bin.cube_pore_sweep_setup:sweep_setup',
              'psm-run-one-atom-sweep = htsohm.bin.one_atom_sweep_run:run_materials',
          ]
      },
)
