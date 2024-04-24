.. _commands/classify_wf:

classify_wf
===========

Classify workflow
-----------------

For arguments and output files, see each of the individual steps:

* :ref:`commands/identify`
* :ref:`commands/align`
* :ref:`commands/classify`

The classify workflow consists of four steps: ``ani_screen``, ``identify``, ``align``, and ``classify``.

The ``ani_screen`` step compares user genomes against a `Mash <https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0997-x>`_ database composed of all GTDB representative genomes,
then verify the best mash hits using `skani <https://www.nature.com/articles/s41592-023-02018-3>`_. User genomes classified with skani are not run through the rest of the pipeline (``identify``, ``align``, ``classify``)
and are reported in the summary file.

The ``identify`` step calls genes using `Prodigal <http://compbio.ornl.gov/prodigal/>`_,
and uses HMM models and the `HMMER <http://hmmer.org/>`_ package to identify the
120 bacterial and 53 archaeal marker genes used for phylogenetic inference
(`Parks et al., 2018 <https://www.ncbi.nlm.nih.gov/pubmed/30148503>`_). Multiple
sequence alignments (MSA) are obtained by aligning marker genes to their respective HMM model.


The ``align`` step concatenates the aligned marker genes and filters the concatenated MSA to approximately 5,000 amino acids.


Finally, the ``classify`` step uses `pplacer <http://matsen.fhcrc.org/pplacer/>`_ to find the maximum-likelihood
placement of each genome in the GTDB-Tk reference tree. GTDB-Tk classifies each genome based on its placement in the reference tree,
its relative evolutionary divergence, and/or average nucleotide identity (ANI) to reference genomes.

Results can be impacted by a lack of marker genes or contamination. We have validated GTDB-Tk on genomes
estimated to be ≥50% complete with ≤10% contamination consistent with community standards for medium or
higher quality single-amplified and metagenome-assembled genomes (`Bowers et al., 2017 <https://www.ncbi.nlm.nih.gov/pubmed/28787424>`_).


The classify workflow can be run as follows:

.. code-block:: bash

    gtdbtk classify_wf --genome_dir <my_genomes> --out_dir <output_dir>

This will process all genomes in the directory <my_genomes> using both bacterial and archaeal marker sets and place the results in <output_dir>. Genomes must be in FASTA format (gzip with the extension .gz is acceptable).
The location of genomes can also be specified using a batch file with the ``--batchfile`` flag. The batch file is a two column file indicating the location of each genome and the desired genome identifier
(i.e., a Newick compatible alphanumeric string). These fields must be separated by a tab.

The workflow supports several optional flags, including:

* ``min_perc_aa``: allows filtering of genomes below a specified percentage of amino acids in the MSA
* ``cpus``: maximum number of CPUs to use

The taxonomic classification of each bacterial and archaeal genome is contained in the
:ref:`[prefix].[domain].summary.tsv <files/summary.tsv>`  output files.

For other flags please consult the command line interface.
 
Arguments
---------

.. argparse::
   :module: gtdbtk.cli
   :func: get_main_parser
   :prog: gtdbtk
   :path: classify_wf
   :nodefaultconst:


Example
-------

Input
^^^^^

.. code-block:: bash

    gtdbtk classify_wf --genome_dir genomes/ --out_dir classify_wf_out --cpus 3



Output
^^^^^^


.. code-block:: text

    [2024-03-25 12:25:02] INFO: GTDB-Tk v2.3.2
    [2024-03-25 12:25:02] INFO: gtdbtk classify_wf --batchfile genomes/5K_batchfile.tsv --out_dir 5K_classify_wf --cpus 64 --mash_db mash_db.msh
    [2024-03-25 12:25:02] INFO: Using GTDB-Tk reference data version r214: /srv/db/gtdbtk/official/release214_skani/release214
    [2024-03-25 12:25:02] INFO: Loading reference genomes.
    [2024-03-25 12:25:03] INFO: Using Mash version 2.2.2
    [2024-03-25 12:25:03] INFO: Creating Mash sketch file: 5K_classify_wf/classify/ani_screen/intermediate_results/mash/gtdbtk.user_query_sketch.msh
    [2024-03-25 12:25:12] INFO: Completed 5,000 genomes in 9.35 seconds (534.79 genomes/second).
    [2024-03-25 12:25:12] INFO: Loading data from existing Mash sketch file: mash_db.msh
    [2024-03-25 12:25:16] INFO: Calculating Mash distances.
    [2024-03-25 12:33:09] INFO: Calculating ANI with skani v0.2.1.
    [2024-03-25 12:34:26] INFO: Completed 40,673 comparisons in 1.26 minutes (32,186.56 comparisons/minute).
    [2024-03-25 12:39:22] INFO: Summary of results saved to: 5K_classify_wf/classify/ani_screen/gtdbtk.bac120.ani_summary.tsv
    [2024-03-25 12:39:22] INFO: 3528 genome(s) have been classified using the ANI pre-screening step.
    [2024-03-25 12:39:22] INFO: Done.
    [2024-03-25 12:39:23] INFO: Identifying markers in 1,472 genomes with 64 threads.
    [2024-03-25 12:39:23] TASK: Running Prodigal V2.6.3 to identify genes.
    [2024-03-25 12:45:33] INFO: Completed 1,472 genomes in 6.17 minutes (238.71 genomes/minute).
    [2024-03-25 12:45:33] TASK: Identifying TIGRFAM protein families.
    [2024-03-25 12:47:22] INFO: Completed 1,472 genomes in 1.82 minutes (809.77 genomes/minute).
    [2024-03-25 12:47:22] TASK: Identifying Pfam protein families.
    [2024-03-25 12:47:29] INFO: Completed 1,472 genomes in 6.77 seconds (217.31 genomes/second).
    [2024-03-25 12:47:29] INFO: Annotations done using HMMER 3.1b2 (February 2015).
    [2024-03-25 12:47:29] TASK: Summarising identified marker genes.
    [2024-03-25 12:48:02] INFO: Completed 1,472 genomes in 32.99 seconds (44.62 genomes/second).
    [2024-03-25 12:48:02] INFO: Done.
    [2024-03-25 12:48:03] INFO: Aligning markers in 1,472 genomes with 64 CPUs.
    [2024-03-25 12:48:03] INFO: Processing 1,472 genomes identified as bacterial.
    [2024-03-25 12:48:15] INFO: Read concatenated alignment for 80,789 GTDB genomes.
    [2024-03-25 12:48:15] TASK: Generating concatenated alignment for each marker.
    [2024-03-25 12:48:23] INFO: Completed 1,472 genomes in 1.38 seconds (1,068.13 genomes/second).
    [2024-03-25 12:48:24] TASK: Aligning 120 identified markers using hmmalign 3.1b2 (February 2015).
    [2024-03-25 12:48:50] INFO: Completed 120 markers in 18.95 seconds (6.33 markers/second).
    [2024-03-25 12:48:51] TASK: Masking columns of bacterial multiple sequence alignment using canonical mask.
    [2024-03-25 12:51:17] INFO: Completed 82,261 sequences in 2.42 minutes (33,935.16 sequences/minute).
    [2024-03-25 12:51:17] INFO: Masked bacterial alignment from 41,084 to 5,035 AAs.
    [2024-03-25 12:51:17] INFO: 0 bacterial user genomes have amino acids in <10.0% of columns in filtered MSA.
    [2024-03-25 12:51:17] INFO: Creating concatenated alignment for 82,261 bacterial GTDB and user genomes.
    [2024-03-25 12:51:42] INFO: Creating concatenated alignment for 1,472 bacterial user genomes.
    [2024-03-25 12:51:42] INFO: Done.
    [2024-03-25 12:51:44] TASK: Placing 1,472 bacterial genomes into backbone reference tree with pplacer using 64 CPUs (be patient).
    [2024-03-25 12:51:44] INFO: pplacer version: v1.1.alpha19-0-g807f6f3
    [2024-03-25 12:56:39] INFO: Calculating RED values based on reference tree.
    [2024-03-25 12:56:45] INFO: 1472 out of 1472 have an class assignments. Those genomes will be reclassified.
    [2024-03-25 12:56:45] TASK: Placing 316 bacterial genomes into class-level reference tree 7 (1/8) with pplacer using 64 CPUs (be patient).
    [2024-03-25 13:03:32] INFO: Calculating RED values based on reference tree.
    [2024-03-25 13:03:34] TASK: Traversing tree to determine classification method.
    [2024-03-25 13:03:34] INFO: Completed 316 genomes in 0.13 seconds (2,474.33 genomes/second).
    [2024-03-25 13:03:35] TASK: Calculating average nucleotide identity using skani (v0.2.1).
    [2024-03-25 13:03:59] INFO: Completed 3,690 comparisons in 24.41 seconds (151.18 comparisons/second).
    [2024-03-25 13:04:26] INFO: 0 genome(s) have been classified using skani and pplacer.
    [2024-03-25 13:04:26] TASK: Placing 284 bacterial genomes into class-level reference tree 3 (2/8) with pplacer using 64 CPUs (be patient).
    [2024-03-25 13:13:22] INFO: Calculating RED values based on reference tree.
    [2024-03-25 13:13:25] TASK: Traversing tree to determine classification method.
    [2024-03-25 13:13:26] INFO: Completed 284 genomes in 0.60 seconds (477.15 genomes/second).
    [2024-03-25 13:13:28] TASK: Calculating average nucleotide identity using skani (v0.2.1).
    [2024-03-25 13:14:12] INFO: Completed 6,618 comparisons in 43.91 seconds (150.72 comparisons/second).
    [2024-03-25 13:14:37] INFO: 0 genome(s) have been classified using skani and pplacer.
    [2024-03-25 13:14:38] TASK: Placing 201 bacterial genomes into class-level reference tree 1 (3/8) with pplacer using 64 CPUs (be patient).
    [2024-03-25 13:22:53] INFO: Calculating RED values based on reference tree.
    [2024-03-25 13:22:56] TASK: Traversing tree to determine classification method.
    [2024-03-25 13:22:56] INFO: Completed 201 genomes in 0.30 seconds (673.19 genomes/second).
    [2024-03-25 13:22:57] TASK: Calculating average nucleotide identity using skani (v0.2.1).
    [2024-03-25 13:23:15] INFO: Completed 2,230 comparisons in 17.76 seconds (125.60 comparisons/second).
    [2024-03-25 13:23:31] INFO: 0 genome(s) have been classified using skani and pplacer.
    [2024-03-25 13:23:31] TASK: Placing 174 bacterial genomes into class-level reference tree 4 (4/8) with pplacer using 64 CPUs (be patient).
    [2024-03-25 13:30:40] INFO: Calculating RED values based on reference tree.
    [2024-03-25 13:30:42] TASK: Traversing tree to determine classification method.
    [2024-03-25 13:30:43] INFO: Completed 174 genomes in 0.70 seconds (249.85 genomes/second).
    [2024-03-25 13:30:44] TASK: Calculating average nucleotide identity using skani (v0.2.1).
    [2024-03-25 13:31:14] INFO: Completed 4,185 comparisons in 29.77 seconds (140.56 comparisons/second).
    [2024-03-25 13:31:29] INFO: 0 genome(s) have been classified using skani and pplacer.
    [2024-03-25 13:31:30] TASK: Placing 170 bacterial genomes into class-level reference tree 6 (5/8) with pplacer using 64 CPUs (be patient).
    [2024-03-25 13:38:23] INFO: Calculating RED values based on reference tree.
    [2024-03-25 13:38:25] TASK: Traversing tree to determine classification method.
    [2024-03-25 13:38:26] INFO: Completed 170 genomes in 0.35 seconds (485.99 genomes/second).
    [2024-03-25 13:38:26] TASK: Calculating average nucleotide identity using skani (v0.2.1).
    [2024-03-25 13:38:44] INFO: Completed 2,192 comparisons in 17.03 seconds (128.74 comparisons/second).
    [2024-03-25 13:38:58] INFO: 0 genome(s) have been classified using skani and pplacer.
    [2024-03-25 13:38:58] TASK: Placing 139 bacterial genomes into class-level reference tree 2 (6/8) with pplacer using 64 CPUs (be patient).
    [2024-03-25 13:47:20] INFO: Calculating RED values based on reference tree.
    [2024-03-25 13:47:23] TASK: Traversing tree to determine classification method.
    [2024-03-25 13:47:24] INFO: Completed 139 genomes in 0.17 seconds (814.69 genomes/second).
    [2024-03-25 13:47:24] TASK: Calculating average nucleotide identity using skani (v0.2.1).
    [2024-03-25 13:47:41] INFO: Completed 1,854 comparisons in 15.93 seconds (116.40 comparisons/second).
    [2024-03-25 13:47:51] INFO: 0 genome(s) have been classified using skani and pplacer.
    [2024-03-25 13:47:52] TASK: Placing 123 bacterial genomes into class-level reference tree 8 (7/8) with pplacer using 64 CPUs (be patient).
    [2024-03-25 13:50:55] INFO: Calculating RED values based on reference tree.
    [2024-03-25 13:50:56] TASK: Traversing tree to determine classification method.
    [2024-03-25 13:50:56] INFO: Completed 123 genomes in 0.02 seconds (5,030.07 genomes/second).
    [2024-03-25 13:50:57] TASK: Calculating average nucleotide identity using skani (v0.2.1).
    [2024-03-25 13:51:06] INFO: Completed 651 comparisons in 8.51 seconds (76.49 comparisons/second).
    [2024-03-25 13:51:17] INFO: 0 genome(s) have been classified using skani and pplacer.
    [2024-03-25 13:51:17] TASK: Placing 65 bacterial genomes into class-level reference tree 5 (8/8) with pplacer using 64 CPUs (be patient).
    [2024-03-25 13:57:05] INFO: Calculating RED values based on reference tree.
    [2024-03-25 13:57:08] TASK: Traversing tree to determine classification method.
    [2024-03-25 13:57:08] INFO: Completed 65 genomes in 0.15 seconds (436.27 genomes/second).
    [2024-03-25 13:57:08] TASK: Calculating average nucleotide identity using skani (v0.2.1).
    [2024-03-25 13:57:18] INFO: Completed 799 comparisons in 9.03 seconds (88.47 comparisons/second).
    [2024-03-25 13:57:24] INFO: 0 genome(s) have been classified using skani and pplacer.
    [2024-03-25 13:57:24] WARNING: 14 of 5000 genomes have a warning (see summary file).
    [2024-03-25 13:57:24] INFO: Note that Tk classification mode is insufficient for publication of new taxonomic designations. New designations should be based on one or more de novo trees, an example of which can be produced by Tk in de novo mode.
    [2024-03-25 13:57:24] INFO: Done.
    [2024-03-25 13:57:24] INFO: Removing intermediate files.
    [2024-03-25 13:57:43] INFO: Intermediate files removed.
    [2024-03-25 13:57:43] INFO: Done.
