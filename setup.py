import os
import re

from setuptools import setup, find_packages


def read_meta():
    """Read each of the keys stored in __init__.py

    Returns
    -------
    dict[str, str]
        A dictionary containing each of the string key value pairs.
    """
    path = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'gtdbtk/__init__.py')
    with open(path) as fh:
        hits = re.findall(r'__(\w+)__ ?= ?["\'](.+)["\']\n', fh.read())
    return {k: v for k, v in hits}


def readme():
    with open('README.md') as f:
        return f.read()


meta = read_meta()
setup(
    author=meta['author'],
    author_email=meta['author_email'],
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
    data_files=[("", ["LICENSE"])],
    description=meta['description'],
    entry_points={
        'console_scripts': [
            'gtdbtk = gtdbtk.__main__:main'
        ]
    },
    install_requires=["dendropy>=4.1.0", 'numpy>=1.9.0', 'tqdm>=4.31.0'],
    license=meta['license'],
    long_description=readme(),
    long_description_content_type='text/markdown',
    maintainer=meta['maintainer'],
    maintainer_email=meta['maintainer_email'],
    name=meta['name'],
    packages=find_packages(),
    package_data={'gtdbtk': ['VERSION', 'tests/data/genomes/*']},
    python_requires=meta['python_requires'],
    url=meta['url'],
    version=meta['version']
)
