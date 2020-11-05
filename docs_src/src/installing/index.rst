.. _installing:

Installing GTDB-Tk
==================

GTDB-Tk is available through multiple sources.

If you are unsure which source to install, Bioconda is generally the easiest.


Sources
-------


.. toctree::
   :maxdepth: 1

   pip
   bioconda
   docker

Alternatively, GTDB-Tk can be run online through `KBase <https://kbase.us/applist/apps/kb_gtdbtk/run_kb_gtdbtk>`_ (third party).


Hardware requirements
---------------------

.. list-table::
   :widths: 10 10 10 30
   :header-rows: 1

   * - Domain
     - Memory
     - Storage
     - Time
   * - Archaea
     - ~8 GB
     - ~27 GB
     - ~1 hour / 1,000 genomes @ 64 CPUs
   * - Bacteria
     - ~150 GB
     - ~27 GB
     - ~1 hour / 1,000 genomes @ 64 CPUs

.. note::
   The amount reported of memory reported can vary depending on the number of pplacer threads.
   See :ref:`faq_pplacer` for more information.


Python libraries
----------------

GTDB-Tk is designed for Python >=3.6 and requires the following libraries, which will be automatically installed:

.. list-table::
   :widths: 10 10 80
   :header-rows: 1

   * - Library
     - Version
     - Reference
   * - `DendroPy <https://dendropy.org/>`_
     - >= 4.1.0
     - Sukumaran, J. and Mark T. Holder. 2010. DendroPy: A Python library for phylogenetic computing. Bioinformatics 26: 1569-1571.
   * - `NumPy <https://numpy.org/>`_
     - >= 1.9.0
     - Harris, C.R., Millman, K.J., van der Walt, S.J. et al. Array programming with NumPy. Nature 585, 357–362 (2020). DOI: `0.1038/s41586-020-2649-2 <https://doi.org/10.1038/s41586-020-2649-2>`_
   * - `tqdm <https://github.com/tqdm/tqdm>`_
     - >= 4.31.0
     - DOI: `10.5281/zenodo.595120 <https://doi.org/10.5281/zenodo.595120>`_


Please cite these libraries if you use GTDB-Tk in your work.


.. _installing#third-party-software:

Third-party software
--------------------

GTDB-Tk makes use of the following 3rd party dependencies and assumes they are on your system path:

.. tip::
   The :ref:`commands/check_install` command will verify that all of the programs are on the path.

.. list-table::
   :widths: 10 10 80
   :header-rows: 1

   * - Software
     - Version
     - Reference
   * - `Prodigal <http://compbio.ornl.gov/prodigal/>`_
     - >= 2.6.2
     - Hyatt D, et al. 2010. `Prodigal: prokaryotic gene recognition and translation initiation site identification <https://www.ncbi.nlm.nih.gov/pubmed/20211023>`_. *BMC Bioinformatics*, 11:119. doi: 10.1186/1471-2105-11-119.
   * - `HMMER <http://hmmer.org/>`_
     - >= 3.1b2
     - Eddy SR. 2011. `Accelerated profile HMM searches <https://www.ncbi.nlm.nih.gov/pubmed/22039361>`_. *PLOS Comp. Biol.*, 7:e1002195.
   * - `pplacer <http://matsen.fhcrc.org/pplacer/>`_
     - >= 1.1
     - Matsen FA, et al. 2010. `pplacer: linear time maximum-likelihood and Bayesian phylogenetic placement of sequences onto a fixed reference tree <https://www.ncbi.nlm.nih.gov/pubmed/21034504>`_. *BMC Bioinformatics*, 11:538.
   * - `FastANI <https://github.com/ParBLiSS/FastANI>`_
     - >= 1.32
     - Jain C, et al. 2019. `High-throughput ANI Analysis of 90K Prokaryotic Genomes Reveals Clear Species Boundaries <https://www.nature.com/articles/s41467-018-07641-9>`_. *Nat. Communications*, doi: 10.1038/s41467-018-07641-9.
   * - `FastTree <http://www.microbesonline.org/fasttree/>`_
     - >= 2.1.9
     - Price MN, et al. 2010. `FastTree 2 - Approximately Maximum-Likelihood Trees for Large Alignments <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2835736/>`_. *PLoS One*, 5, e9490.
   * - `Mash <https://github.com/marbl/Mash>`_
     - >= 2.2
     - Ondov BD, et al. 2016. `Mash: fast genome and metagenome distance estimation using MinHash <https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0997-x>`_. *Genome Biol* 17, 132. doi: doi: 10.1186/s13059-016-0997-x.


Please cite these tools if you use GTDB-Tk in your work.


.. _installing#gtdbtk-reference-data:

GTDB-Tk reference data
----------------------

GTDB-Tk requires ~27G of external data that needs to be downloaded and unarchived:

.. code-block:: bash

    wget https://data.ace.uq.edu.au/public/gtdb/data/releases/release95/95.0/auxillary_files/gtdbtk_r95_data.tar.gz
    tar xvzf gtdbtk_r95_data.tar.gz


Note that different versions of the GTDB release data may not run on all versions of GTDB-Tk, below are all supported versions:


.. list-table::
   :widths: 10 10 10
   :header-rows: 1

   * - GTDB Release
     - Minimum version
     - Maximum version
   * - R95
     - 1.3.0
     - N/A
   * - R89
     - 0.3.0
     - 0.1.2
   * - R86.2
     - 0.2.1
     - 0.2.2
   * - R86
     - 0.1.0
     - 0.1.6
   * - R83
     - 0.0.6
     - 0.0.7


Reference data for prior releases of GTDB-Tk are available at: https://data.ace.uq.edu.au/public/gtdbtk
