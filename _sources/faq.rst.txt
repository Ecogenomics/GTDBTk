.. _faq:

FAQ
===

Taxonomy FAQ
------------

Why is there a discrepancy in the naming system between GTDB-Tk and NCBI or Silva taxonomic names?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

GTDB-Tk uses the GTDB taxonomy (`https://gtdb.ecogenomic.org/ <https://gtdb.ecogenomic.org/>`_).
This taxonomy is similar, but not identical to NCBI and Silva.
In many cases the GTDB taxonomy more strictly follows the nomenclatural rules for rank suffixes which is why there is Nitrospirota instead of Nitrospirae.


Can you combine the bacterial and archaeal trees into a single tree?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


The bacterial and archaeal trees are inferred from different marker genes. Currently, the correct rootings of these trees remain an open area of research.
GTDB-Tk does not provide a tool to merge the trees but It is possible to artificially combine them by manipulating the Newick files.
One solution would be to use (`DendroPy <https://dendropy.org/>`_); a Python library used for phylogenetic computing.


GTDB-Tk FAQ
------------

.. _faq_pplacer:

GTDB-Tk reaches the memory limit / pplacer crashes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The host may report that GTDB-Tk has exceeded the memory requirements due to how ``pplacer`` is implemented.
Briefly, this is only the reported value and is not true for how much memory is actually in use.

When pplacer runs, it goes through several steps notable detailed below:

#. pplacer requests an uninitialised chunk of memory from the host (say, 150 GB).

   * VIRT = 150 GB ``[virtual memory memory, i.e. the amount of memory mapped]``

   * RES ~= 0 GB ``[amount of memory physically in use]``

#. pplacer starts caching information to that chunk of memory

   * VIRT = 150 GB

   * RES -> 150 GB  (slowly increases to 150 GB as caching information is written)

   .. note::
      If pplacer crashes here, you likely don't have enough free memory on the server.

#. pplacer starts placing genomes into the tree

   #. the main pplacer thread forks for each ``--pplacer_cpus`` or ``--cpus`` (if not specified).

      * Unix uses the `copy-on-write <https://en.wikipedia.org/wiki/Copy-on-write>`_ method for each child.

      .. note::
         This means that a new thread will appear using the same amount of memory as the parent (150 GB).
         Due to how copy-on-write is implemented, this is the same memory space and is not using any additional physical memory.

         Therefore, the host will report pplacer using a total of ``PARENT_MEMORY * (N_CHILDREN + 1)`` GB.

   #. each worker will only read from the memory space and exit once the queue of query genomes is depleted.


For example, running GTDB-Tk with on the bacterial tree (requires 150 GB of memory) with 1 CPU will require 150 GB of physical
memory, but the host will report 300 GB of memory in use.

Using the ``--scratch_dir`` parameter and ``--pplacer_cpus 1`` may help.


How is GTDB-Tk validating species assignments using average nucleotide identity?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

GTDB-Tk uses `FastANI <https://github.com/ParBLiSS/FastANI>`_ to estimate the ANI between genomes.
We recommend you have FastANI >= 1.32 as this version introduces a fix that makes the results deterministic.
A query genome is only classified as belonging to the same species as a reference genome if the ANI between the
genomes is within the species ANI circumscription radius (typically, 95%) and the alignment fraction (AF) is >=0.5.
In some circumstances, the phylogenetic placement of a query genome may not support the species assignment.
GTDB r207 strictly uses ANI to circumscribe species and GTDB-Tk follows this methodology.
The species-specific ANI circumscription radii are available from the `GTDB <https://gtdb.ecogenomic.org/>`_ website.


Why is FastANI using more threads than allocated?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you are using FastANI version 1.33 then you may run into an issue where FastANI will use more threads than you allocate.
This can be problematic if running GTDB-Tk on a HPC where you have a limited number of threads available.

This issue has been `reported to the FastANI developers (#101) <https://github.com/ParBLiSS/FastANI/issues/101>`_.

Depending on how you installed GTDB-Tk there are different ways to downgrade FastANI to version 1.32.

**Manual:**

Simplify download and install the FastANI binary from `here <https://github.com/ParBLiSS/FastANI/releases/tag/v1.32>`_.

**Conda:**

From GTDB-Tk v2.0.0 the conda environment will automatically have FastANI v1.3 installed. Otherwise run:

``conda install -c bioconda fastani==1.32``

**Docker:**

From GTDB-Tk v2.2.2 the Docker container will automatically have FastANI v1.32 installed. Otherwise, manually
build the container from the `Dockerfile <https://github.com/Ecogenomics/GTDBTk/blob/master/Dockerfile>`_, making
sure to specify FastANI v1.32.

What is the difference between the mutually exclusive options ``--mash_db`` and ``--skip_ani_screen``?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

| Starting with GTDB-Tk v2.2+, the ``classify_wf`` and ``classify`` function require an extra parameter to run: ``--mash_db`` or ``--skip_ani_screen``.
| With this new version of Tk, The first stage of ``classify`` pipelines (``classify_wf`` and ``classify``) is to compare all user genomes to all reference genomes and annotate them, if possible, based on ANI matches.
| Using the ``--mash_db`` option will indicate to GTDB-Tk the path of the sketched Mash database require for ANI screening.
| If no database are available ( i.e. this is the first time running classify ), the ``--mash_db`` option will sketch a new Mash database that can be used for subsequent calls.
| The ``--skip_ani_screen`` option will skip the pre-screening step and classify all genomes similar to previous versions of GTDB-Tk.