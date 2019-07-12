#!/usr/bin/env python

import os

from setuptools import setup


def version():
    setup_dir = os.path.dirname(os.path.realpath(__file__))
    with open(os.path.join(setup_dir, 'gtdbtk', 'VERSION'), 'r') as f:
        return f.readline().strip()


setup(
    name='gtdbtk',
    version=version(),
    author='Pierre-Alain Chaumeil and Donovan Parks',
    author_email='p.chaumeil@uq.edu.au',
    maintainer='Pierre-Alain Chaumeil and Donovan Parks',
    maintainer_email='donovan.parks@gmail.com',
    packages=['gtdbtk', 'gtdbtk.biolib_lite', 'gtdbtk.config', 'gtdbtk.external', 'gtdbtk.tests',
              'gtdbtk.external.pypfam', 'gtdbtk.external.pypfam.HMM', 'gtdbtk.external.pypfam.Scan'],
    scripts=['bin/gtdbtk'],
    package_data={'gtdbtk': ['VERSION',
                             'MANIFEST.in', 'tests/data/genomes/*']},
    url='http://pypi.python.org/pypi/gtdbtk/',
    license='GPL3',
    description='A toolkit for assigning objective taxonomic classifications to bacterial and archaeal genomes.',
    classifiers=[
        'Development Status :: 4 - Beta',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Natural Language :: English',
        'Programming Language :: Python :: 2.7',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    install_requires=["dendropy>=4.1.0", 'numpy', 'future'],
)
