.. _installing:

Installing GTDB-Tk
==================

GTDB-Tk is available through multiple sources.

If you are unsure which source to install, Bioconda is generally the easiest.


Sources
-------

Local
^^^^^

.. toctree::
   :maxdepth: 1

   pip
   bioconda
   docker


Online
^^^^^^

.. toctree::
   :maxdepth: 1

   kbase


Hardware requirements
---------------------

* ~100Gb of memory to run
* ~27Gb of storage
* ~1 hour per 1,000 genomes when using 64 CPUs


Python libraries
----------------

GTDB-Tk is designed for Python >=3.6 and requires the following libraries, which will be automatically installed:

* `dendropy <http://dendropy.org/>`_  >=4.1.0: a Python library for phylogenetic computing.
* `NumPy <https://numpy.org/>`_ >=1.9.0: scientific computing with Python.
* `tqdm <https://github.com/tqdm/tqdm>`_: A Fast, Extensible Progress Bar for Python and CLI https://tqdm.github.io


Third-party software
--------------------

GTDB-Tk makes use of the following 3rd party dependencies and assumes they are on your system path:

* `Prodigal <http://compbio.ornl.gov/prodigal/>`_ >= 2.6.2: Hyatt D, et al. 2010. `Prodigal: prokaryotic gene recognition and translation initiation site identification <https://www.ncbi.nlm.nih.gov/pubmed/20211023>`_.
  *BMC Bioinformatics*, 11:119. doi: 10.1186/1471-2105-11-119.

* `HMMER <http://hmmer.org/>`_ >= 3.1: Eddy SR. 2011. `Accelerated profile HMM searches <https://www.ncbi.nlm.nih.gov/pubmed/22039361>`_.
  *PLOS Comp. Biol.*, 7:e1002195.

* `pplacer <http://matsen.fhcrc.org/pplacer/>`_ >= 1.1: Matsen FA, et al. 2010. `pplacer: linear time maximum-likelihood and Bayesian phylogenetic placement of sequences onto a fixed reference tree <https://www.ncbi.nlm.nih.gov/pubmed/21034504>`_.
  *BMC Bioinformatics*, 11:538.

* `FastANI <https://github.com/ParBLiSS/FastANI>`_ >= 1.2: Jain C, et al. 2019. `High-throughput ANI Analysis of 90K Prokaryotic Genomes Reveals Clear Species Boundaries <https://www.nature.com/articles/s41467-018-07641-9>`_.
  *Nat. Communications*, doi: 10.1038/s41467-018-07641-9.

* `FastTree <http://www.microbesonline.org/fasttree/>`_ >= 2.1.9: Price MN, et al. 2010. `FastTree 2 - Approximately Maximum-Likelihood Trees for Large Alignments <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2835736/>`_.
  *PLoS One*, 5, e9490.

* `Mash <https://github.com/marbl/Mash>`_ >= 2.2: Ondov BD, et al. 2016. `Mash: fast genome and metagenome distance estimation using MinHash <https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0997-x>`_.
  *Genome Biol* 17, 132. doi: doi: 10.1186/s13059-016-0997-x.

Please cite these tools if you use GTDB-Tk in your work.

GTDB-Tk reference data
----------------------

GTDB-Tk requires ~27G of external data that need to be downloaded and unarchived:

.. code-block:: bash

    wget https://data.ace.uq.edu.au/public/gtdb/data/releases/release95/95.0/auxillary_files/gtdbtk_r95_data.tar.gz
    tar xvzf gtdbtk_r95_data.tar.gz


Reference data for prior releases of GTDB-Tk are available at: https://data.ace.uq.edu.au/public/gtdbtk

