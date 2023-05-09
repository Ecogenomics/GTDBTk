.. _installing:

Installing GTDB-Tk
==================

GTDB-Tk is available through multiple sources, you only need to choose one.
If you are unsure which one to choose, Bioconda is generally the easiest.


Sources
-------


.. toctree::
   :maxdepth: 1

   bioconda
   pip
   docker

Alternatively, GTDB-Tk can be run online through `KBase <https://kbase.us/applist/apps/kb_gtdbtk/run_kb_gtdbtk>`_ (third party). Note that the version may not be the most recent release.


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
     - ~45 GB
     - ~85 GB
     - ~1 hour / 1,000 genomes @ 64 CPUs
   * - Bacteria
     - ~65GB (410 GB when using --full_tree)
     - ~85 GB
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
     - >= 4.35.0
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

GTDB-Tk requires ~84G of external data that needs to be downloaded and unarchived:

.. code-block:: bash

    wget https://data.gtdb.ecogenomic.org/releases/latest/auxillary_files/gtdbtk_data.tar.gz
    wget https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/auxillary_files/gtdbtk_data.tar.gz  (or, mirror)
    tar xvzf gtdbtk_data.tar.gz


.. note:: Note that different versions of the GTDB release data may not run on all versions of GTDB-Tk, check the supported versions!


.. list-table::
   :widths: 10 10 10 20
   :header-rows: 1

   * - GTDB Release
     - Minimum version
     - Maximum version
     - MD5
   * - `R214 <https://data.gtdb.ecogenomic.org/releases/release214/214.0/auxillary_files/gtdbtk_r214_data.tar.gz>`_
     - 2.1.0
     - Current
     - 630745840850c532546996b22da14c27
   * - `R207_v2 <https://data.gtdb.ecogenomic.org/releases/release207/207.0/auxillary_files/gtdbtk_r207_v2_data.tar.gz>`_
     - 2.1.0
     - Current
     - df468d63265e8096d8ca01244cb95f30
   * - `R207 <https://data.gtdb.ecogenomic.org/releases/release207/207.0/auxillary_files/gtdbtk_r207_data.tar.gz>`_
     - 2.0.0
     - 2.0.0
     - b04c55104b491f84e053a9011b36164a
   * - `R202 <https://data.gtdb.ecogenomic.org/releases/release202/202.0/auxillary_files/gtdbtk_r202_data.tar.gz>`_
     - 1.5.0
     - 1.7.0
     - 4986526c2b935fd4dcc2e604c0322517
   * - `R95 <https://data.gtdb.ecogenomic.org/releases/release95/95.0/auxillary_files/gtdbtk_r95_data.tar.gz>`_
     - 1.3.0
     - 1.4.2
     - 06924c63f4b555ac6fd1525b09901186
   * - `R89 <https://data.gtdb.ecogenomic.org/releases/release89/89.0/gtdbtk_r89_data.tar.gz>`_
     - 0.3.0
     - 0.1.2
     - 82966ef36086237d7230955e2bfff759
   * - `R86.2 <https://data.gtdb.ecogenomic.org/releases/release86/86.2/gtdbtk.r86_v2_data.tar.gz>`_
     - 0.2.1
     - 0.2.2
     - f71408d69fa2a289f2cdc734b7a58a02
   * - `R86 <https://data.gtdb.ecogenomic.org/releases/release86/86.0/gtdbtk_r86_data.tar.gz>`_
     - 0.1.0
     - 0.1.6
     - d019b3541746c3673181f24e666594ba
   * - `R83 <https://data.gtdb.ecogenomic.org/releases/release83/83.0/gtdbtk_r83_data.tar.gz>`_
     - 0.0.6
     - 0.0.7
     - 9cf523761da843b5787f591f6c5a80de
