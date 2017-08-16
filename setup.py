from setuptools import setup

import os

def version():
    setupDir = os.path.dirname(os.path.realpath(__file__))
    versionFile = open(os.path.join(setupDir, 'src', 'VERSION'))
    return versionFile.read().strip()

setup(
    name='gtdbtk',
    version=version(),
    author='Pierre-Alain Chaumeil and Donovan Parks',
    author_email='p.chaumeil@uq.edu.au',
    packages=['gtdbtk', 'gtdbtk.config', 'gtdbtk.external'],
    scripts=['bin/gtdbtk'],
    package_data={'gtdbtk': ['VERSION']},
    url='http://pypi.python.org/pypi/gtdbtk/',
    license='GPL3',
    description='A toolkit for assigning objective taxonomic classifications to bacterial and archaeal genomes.',
    long_description=open('README.md').read(),
    install_requires=[],
)
