#!/usr/bin/env python

import os

from setuptools import setup


def version():
    setup_dir = os.path.dirname(os.path.realpath(__file__))
    with open(os.path.join(setup_dir, 'gtdbtk', 'VERSION'), 'r') as f:
        return f.readline().strip()


def readme():
    with open('README.md') as f:
        return f.read()


setup(
    name='gtdbtk',
    python_requires='>=3.6',
    version=version(),
    author='Pierre-Alain Chaumeil and Donovan Parks',
    author_email='p.chaumeil@uq.edu.au',
    maintainer='Pierre-Alain Chaumeil, Aaron Mussig, and Donovan Parks',
    maintainer_email='donovan.parks@gmail.com',
    packages=['gtdbtk', 'gtdbtk.biolib_lite', 'gtdbtk.config', 'gtdbtk.external', 'gtdbtk.tests',
              'gtdbtk.external.pypfam', 'gtdbtk.external.pypfam.HMM', 'gtdbtk.external.pypfam.Scan',
              'gtdbtk.io', 'gtdbtk.io.marker', 'gtdbtk.io.prodigal'],
    scripts=['bin/gtdbtk'],
    package_data={'gtdbtk': ['VERSION',
                             'MANIFEST.in', 'tests/data/genomes/*']},
    url='https://github.com/Ecogenomics/GTDBTk',
    license='GPL3',
    description='A toolkit for assigning objective taxonomic classifications to bacterial and archaeal genomes.',
    long_description=readme(),
    long_description_content_type='text/markdown',
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    install_requires=["dendropy>=4.1.0", 'numpy'],
    data_files=[("", ["LICENSE"])]
)
