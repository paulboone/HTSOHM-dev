from setuptools import setup

setup(
    name = 'HTSOHM',
    version = '0.6.0',
    description = 'High-throughput Screening of Pseudo Materials',
    author = 'Paul Boone / A Reino Kaija',
    author_email = 'paulboone@pitt.edu',
    url = 'https://github.com/paulboone/htsohm',
    packages = ['htsohm'],
    install_requires=[
        'numpy ~= 1.17',
        'scipy ~= 1.3',
        'sqlalchemy >= 1',
        'click >= 6',
        'pandas ~=0.25',
        'pyyaml >= 3',
    ],
    extras_require={
        'plots': [
            'matplotlib ~= 3.1',
        ],
    },
    tests_requires=[
        "pytest ~= 6.0"
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
              'psm-material-csv = htsohm.bin.output_csv:output_material_csv',
              'psm-csv-add-bin = htsohm.bin.output_csv:csv_add_bin',
              'psm-dof-analysis = htsohm.bin.dof_analysis:dof_analysis',
          ]
      },
)
