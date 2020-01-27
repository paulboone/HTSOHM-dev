#!/usr/bin/env python

from distutils.core import setup
setup(
    name = 'HTSOHM',
    version = '0.6.0',
    description = 'High-throughput Screening of Pseudo Materials',
    author = 'Paul Boone / A Reino Kaija',
    author_email = 'paulboone@pitt.edu',
    url = 'https://github.com/paulboone/htsohm',
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
              'psm-graph-max-prop-by-num-materials = htsohm.bin.graph_max_prop_by_num_materials:bin_graph',
              'psm-graph-prop-by-num-materials = htsohm.bin.graph_prop_by_num_materials:bin_graph',
              'psm-output-config-files = htsohm.bin.config_files:output_config_files',
              'psm-csv = htsohm.bin.output_csv:output_csv',
              'psm-atoms-csv = htsohm.bin.output_csv:output_atom_sites_csv',
              'psm-csv-add-bin = htsohm.bin.output_csv:csv_add_bin',
              'psm-dof-analysis = htsohm.bin.dof_analysis:dof_analysis',
              'psm-setup-one-atom-sweep = htsohm.bin.one_atom_sweep_setup:sweep_setup',
              'psm-setup-cube-pore-sweep = htsohm.bin.cube_pore_sweep_setup:sweep_setup',
              'psm-run-one-atom-sweep = htsohm.bin.one_atom_sweep_run:run_materials',
          ]
      },
)
