.. _commands/de_novo_wf:

de_novo_wf
==========

For arguments and output files, see each of the individual steps:

* :ref:`commands/identify`
* :ref:`commands/align`
* :ref:`commands/infer`
* :ref:`commands/root`
* :ref:`commands/decorate`


The *de novo* workflow infers a new bacterial or archaeal tree containing the user-supplied genomes and, by default,
the GTDB reference genomes. The classify workflow is recommended for obtaining taxonomic classifications; this workflow
is only recommended if a *de novo* domain-specific tree is desired. **The taxonomic assignments should be taken as a
guide, not as final classifications.** In particular, no effort is made to resolve the taxonomic assignment of lineages
composed exclusively of user-submitted genomes.

This workflow consists of five steps: ``identify``, ``align``, ``infer``, ``root`` and ``decorate``.

* ``identify`` calls genes with Prodigal and identifies the marker genes, as in the classify workflow.
  No ANI screen is performed.
* ``align`` builds the multiple sequence alignment of the user genomes and, unless ``--skip_gtdb_refs`` is used, the
  GTDB reference genomes (optionally restricted with ``--taxa_filter``). Columns are selected with the canonical mask,
  or with ``--custom_msa_filters`` and its parameters (``--cols_per_gene``, ``--min_consensus``,
  ``--max_consensus``, ``--min_perc_taxa``, ``--rnd_seed``).
* ``infer`` builds the tree with `FastTree <http://www.microbesonline.org/fasttree/>`_ (FastTreeMP when
  ``--cpus`` > 1) using the ``--prot_model`` substitution model (default: WAG). Branch lengths are rescaled under the
  Gamma20 model only when ``--gamma`` is used, and local support values are computed unless ``--no_support`` is used.
* ``root`` roots the tree on the ``--outgroup_taxon``.
* ``decorate`` decorates the rooted tree with the GTDB taxonomy.

The *de novo* workflow can be run as follows:

.. code-block:: bash

    gtdbtk de_novo_wf --genome_dir <my_genomes> --<bacteria|archaea> --outgroup_taxon <outgroup> --out_dir <output_dir>


This will process all genomes in ``<my_genomes>`` using the specified marker set (``--bacteria`` or ``--archaea``) and
place the results in ``<output_dir>``. Only genomes previously identified as bacterial (archaeal) should be included
when using the bacterial (archaeal) marker set. The tree is rooted on the ``<outgroup>`` taxon (typically a phylum in
the domain-specific tree) as required for correct decoration of the tree. In general, we suggest the resulting tree be
treated as unrooted when interpreting results. As in the classify workflow, genomes can also be specified with a batch
file (``--batchfile``), and gzipped FASTA files can be used with ``--extension gz``.

The workflow supports several optional flags, including:

* ``--cpus``: maximum number of CPUs to use.
* ``--min_perc_aa``: exclude genomes that do not have at least this percentage of amino acids in the MSA
  (default: 10).
* ``--taxa_filter``: restrict the GTDB reference genomes to the given taxa (comma separated, e.g. ``p__Bacillota``).
* ``--skip_gtdb_refs``: do not include GTDB reference genomes in the MSA. Requires ``--custom_taxonomy_file``, which
  must contain the genomes of the outgroup.
* ``--custom_taxonomy_file``: taxonomy of user genomes (see below), used for rooting and decoration.
* ``--prot_model``: protein substitution model for tree inference (``JTT``, ``WAG`` or ``LG``; default: ``WAG``).
* ``--gamma``: rescale branch lengths to optimize the Gamma20 likelihood.
* ``--no_support``: do not compute local support values (Shimodaira-Hasegawa test).
* ``--keep_intermediates``: keep intermediate files in the output directory.

For all flags, see the arguments below or the command line interface.


Arguments
---------

.. argparse::
   :module: gtdbtk.cli
   :func: get_main_parser
   :prog: gtdbtk
   :path: de_novo_wf
   :nodefaultconst:



Example
-------

Input
^^^^^

.. code-block:: bash

    gtdbtk de_novo_wf --genome_dir genomes/ --outgroup_taxon p__Undinarchaeota --archaea --out_dir de_novo_wf --cpus 3

    gtdbtk de_novo_wf --genome_dir genomes/ --outgroup_taxon p__Chloroflexota --bacteria --taxa_filter p__Bacillota,p__Chloroflexota --out_dir de_novo_output

    # Skip GTDB reference genomes (requires --custom_taxonomy_file for the outgroup)
    gtdbtk de_novo_wf --genome_dir genomes/ --outgroup_taxon p__Customphylum --bacteria --skip_gtdb_refs --custom_taxonomy_file custom_taxonomy.tsv --out_dir de_novo_output

    # Use a subset of GTDB reference genomes (p__Bacillota) and root on a custom phylum (p__Customphylum)
    gtdbtk de_novo_wf --genome_dir genomes/ --taxa_filter p__Bacillota --outgroup_taxon p__Customphylum --bacteria --custom_taxonomy_file custom_taxonomy.tsv --out_dir de_novo_output

Custom Taxonomy Format
^^^^^^^^^^^^^^^^^^^^^^
The custom taxonomy file is a tab-separated file with the first column listing user genomes (i.e. the FASTA file name
without the extension) and the second column listing the standardized 7-rank taxonomy.

.. code-block:: text

    # For genome_1.fna, genome_2.fna and genome_3.fna
    genome_1	d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Salmonella;s__Salmonella enterica
    genome_2	d__Bacteria;p__Actinomycetota;c__Actinomycetes;o__Mycobacteriales;f__Mycobacteriaceae;g__Mycobacterium;s__Mycobacterium tuberculosis
    genome_3	d__Bacteria;p__Bacillota;c__Bacilli;o__Lactobacillales;f__Streptococcaceae;g__Streptococcus;s__Streptococcus pyogenes
