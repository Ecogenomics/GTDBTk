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

Why do I see discrepancies in classification for closely related genomes?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Discrepancies in taxonomic assignments can occur when working with closely related genomes. GTDB-Tk uses both the **Relative Evolutionary Divergence (RED)** value and the **placement of the genome in the reference tree** to determine the best taxonomic classification.
In most cases taxonomic assignments are robust, but it is possible for highly similar genomes  to have sufficiently different protein sequences that their placements in the reference tree will vary slightly. This can lead to:

1. placement on different but closely related branches, or
2. placement on the same branch but at different depths.

This uncertainty in the placement of even closely related genomes can result in genomes being erroneously assigned to different taxa (e.g. closely related sister families).

**What should I do?**

If you run into this situation, here are a few strategies you can use:

1. Dereplicate your genomes.
    If you have a cluster of closely related genomes representing a single species, consider dereplicating them to pick a single representative based on genome quality metrics (for example, using a tool like `dRep <https://drep.readthedocs.io/en/latest/index.html>`_ or `Galah <https://github.com/wwood/galah>`_). The taxonomic assignment for this species representative genomes can then be propagated to the other genomes in the cluster.
2. Build a de novo tree
    Another option is to construct a de novo tree including your closely related genomes plus a suitable outgroup. For instance, if you have three genomes classified as:

    * *c__Atribacteria;o__Atribacterales;f__Atribacteraceae*
    * *c__Atribacteria;o__Atribacterales;f__Thermatribacteraceae*
    * *c__Atribacteria;o__Atribacterales;f__*

you could generate a tree for the class Atribacteria (with one outgroup class) and examine whether your genomes cluster together despite their different family-level assignments.

Can you combine the bacterial and archaeal trees into a single tree?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


The bacterial and archaeal trees are inferred from different marker genes. Currently, the correct rootings of these trees remain an open area of research.
GTDB-Tk does not provide a tool to merge the trees but It is possible to artificially combine them by manipulating the Newick files.
One solution would be to use (`DendroPy <https://dendropy.org/>`_); a Python library used for phylogenetic computing.


GTDB-Tk FAQ
------------

.. _faq_memory:

GTDB-Tk reaches the memory limit during placement
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Placing genomes in the reference tree with gPlace is the step that needs the most memory. By default
gPlace places untiled: it holds the likelihood vectors of the whole reference tree in memory, which is
the fastest way to place but needs more than 200 GB for the bacterial tree, more than GTDB-Tk 2.8
needed (~140 GB). gPlace warns in the GTDB-Tk log when its
estimated peak memory is above the memory available to it (the machine's memory, or a SLURM or cgroup limit).

Use ``--gplace_max_memory`` to set the memory gPlace may use, for example ``--gplace_max_memory 64G``
(``K``, ``M``, ``G`` and ``T`` are binary units, as in SLURM; a plain number is GiB). If placement does
not fit within that limit, gPlace places the tree tile by tile: the results are identical, the run is
slower. Leave some memory for the rest of GTDB-Tk, e.g. set it a few GB below your job's memory request.

``--scratch_dir`` from earlier versions has no effect with gPlace, and ``--pplacer_cpus`` is read as
``--gplace_cpus``; both will be removed in the next release.

See the `gPlace documentation <https://cmc-aau.github.io/gPlace/advanced/memory.html>`_ for details.


How is GTDB-Tk validating species assignments using average nucleotide identity?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

GTDB-Tk uses `skani <https://github.com/bluenote-1577/skani>`_ ( it was using fastANI until v2.3.2) to estimate the ANI between genomes.
A query genome is only classified as belonging to the same species as a reference genome if the ANI between the
genomes is within the species ANI circumscription radius (typically, 95%) and the alignment fraction (AF) is >=0.5.
In some circumstances, the phylogenetic placement of a query genome may not support the species assignment.
GTDB r207+ strictly uses ANI to circumscribe species and GTDB-Tk follows this methodology.
The species-specific ANI circumscription radii are available from the `GTDB <https://gtdb.ecogenomic.org/>`_ website.


What is the difference between the mutually exclusive options ``--mash_db`` and ``--skip_ani_screen``?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

| Starting with GTDB-Tk v2.2+, the ``classify_wf`` and ``classify`` function require an extra parameter to run: ``--mash_db`` or ``--skip_ani_screen``.
| With this new version of Tk, The first stage of ``classify`` pipelines (``classify_wf`` and ``classify``) is to compare all user genomes to all reference genomes and annotate them, if possible, based on ANI matches.
| Using the ``--mash_db`` option will indicate to GTDB-Tk the path of the sketched Mash database require for ANI screening.
| If no database are available ( i.e. this is the first time running classify ), the ``--mash_db`` option will sketch a new Mash database that can be used for subsequent calls.
| The ``--skip_ani_screen`` option will skip the pre-screening step and classify all genomes similar to previous versions of GTDB-Tk.

Deprecated FAQ
---------------

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