#!/usr/bin/env python

from setuptools import setup

import os


def version():
    setupDir = os.path.dirname(os.path.realpath(__file__))
    versionFile = open(os.path.join(setupDir, 'gtdbtk', 'VERSION'))
    return versionFile.readline().strip()

setup(
    name='gtdbtk',
    version=version(),
    author='Pierre-Alain Chaumeil and Donovan Parks',
    author_email='p.chaumeil@uq.edu.au',
    maintainer='Pierre-Alain Chaumeil and Donovan Parks',
    maintainer_email='donovan.parks@gmail.com',
    packages=['gtdbtk','gtdbtk.config','gtdbtk.external','gtdbtk.external.Bio.Pfam.Active_site','gtdbtk.external.Bio.Pfam.HMM','gtdbtk.external.Bio.Pfam.Scan'],
    scripts=['bin/gtdbtk'],
    package_data={'gtdbtk': ['VERSION','MANIFEST.in'],
                  'gtdbtk.external': ['pfam_search.pl'],
                  'gtdbtk.external.Bio.Pfam.Active_site':['as_search.pm'],
                  'gtdbtk.external.Bio.Pfam.HMM':['HMM.pm','HMMIO.pm','HMMMatch.pm','HMMResults.pm','HMMResultsIO.pm','HMMSequence.pm','HMMUnit.pm'],
                  'gtdbtk.external.Bio.Pfam.Scan':['PfamScan.pm','Seq.pm'],
                  },
    
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
    install_requires=["biolib >= 0.0.43",
                        "jinja2>=2.7.3",
                        "mpld3>=0.2",
                        "dendropy>=4.1.0"],
)
