#!/usr/bin/env python

from distutils.core import setup
setup(
    name = 'HTSOHM',
    version = '0.4.0',
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
              'dps = htsohm.dps:cmdline'
          ]
      },
)
