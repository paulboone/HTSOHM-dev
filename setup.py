#!/usr/bin/env python

from distutils.core import setup
setup(
    name = 'HTSOHM',
    version = '0.5.0',
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
              'dps = htsohm.dps:cmdline',
              'pm-graph-bins = htsohm.bin.bin_graph:bin_graph',
              'pm-output-config-files = htsohm.bin.config_files:output_config_files',
              'pm-dof-analysis = htsohm.bin.dof_analysis:dof_analysis',
              'pm-setup-one-atom-sweep = htsohm.bin.one_atom_sweep_setup:sweep_materials',
              'pm-run-one-atom-sweep = htsohm.bin.one_atom_sweep_run:run_materials',

          ]
      },
)
