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
then verify the best mash hits using `FastANI <https://www.nature.com/articles/s41467-018-07641-9>`_. User genomes classified with FastANI are not run through the rest of the pipeline (``identify``, ``align``, ``classify``)
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

    [2023-02-22 16:10:50] INFO: GTDB-Tk v2.2.3
    [2023-02-22 16:10:50] INFO: gtdbtk classify_wf --batchfile 3lines_batchfile.tsv --out_dir classify_wf_outdir_test --keep_intermediates --cpus 20 --mash_db mash_sketch/cli/mash_db.msh
    [2023-02-22 16:10:50] INFO: Using GTDB-Tk reference data version r207: /srv/projects/gtdbtk/test_new_features/release207_v2/
    [2023-02-22 16:10:50] INFO: Loading reference genomes.
    [2023-02-22 16:10:51] INFO: Using Mash version 2.3
    [2023-02-22 16:10:51] INFO: Loading data from existing Mash sketch file: classify_wf_outdir_test/classify/ani_screen/intermediate_results/mash/gtdbtk.user_query_sketch.msh
    [2023-02-22 16:10:51] INFO: Creating Mash sketch file: mash_sketch/cli/mash_db.msh
    [2023-02-22 16:10:51] INFO: Calculating RED values based on reference tree.
    [2023-02-22 16:10:54] TASK: Traversing tree to determine classification method.
    [2023-02-22 16:10:54] INFO: Completed 1 genome in 0.00 seconds (2,335.36 genomes/second).
    [2023-02-22 16:10:54] TASK: Calculating average nucleotide identity using FastANI (v1.3).
    [2023-02-22 16:10:57] INFO: Completed 34 comparisons in 2.27 seconds (14.95 comparisons/second).
    [2023-02-22 16:10:57] INFO: 0 genome(s) have been classified using FastANI and pplacer.
    [2023-02-22 16:10:57] TASK: Placing 1 bacterial genomes into class-level reference tree 5 (2/2) with pplacer using 20 CPUs (be patient).
    [2023-02-22 16:14:29] INFO: Calculating RED values based on reference tree.
    [2023-02-22 16:14:31] TASK: Traversing tree to determine classification method.
    [2023-02-22 16:14:31] INFO: Completed 1 genome in 0.06 seconds (16.77 genomes/second).
    [2023-02-22 16:14:31] INFO: 0 genome(s) have been classified using FastANI and pplacer.
    [2023-02-22 16:14:31] WARNING: 1 of 3 genome has a warning (see summary file).
    [2023-02-22 16:14:31] INFO: Note that Tk classification mode is insufficient for publication of new taxonomic designations. New designations should be based on one or more de novo trees, an example of which can be produced by Tk in de novo mode.
    [2023-02-22 16:14:31] INFO: Done.
    [2023-02-22 16:20:06] INFO: Completed 65,703 genomes in 9.25 minutes (7,103.32 genomes/minute).
    [2023-02-22 16:20:06] INFO: Calculating Mash distances.
    [2023-02-22 16:20:10] INFO: Calculating ANI with FastANI v1.3.
    [2023-02-22 16:20:11] INFO: Completed 12 comparisons in 0.63 seconds (18.90 comparisons/second).
    [2023-02-22 16:20:11] INFO: Summary of results saved to: classify_wf_outdir_test_mash/classify/ani_screen/gtdbtk.bac120.ani_summary.tsv
    [2023-02-22 16:20:11] INFO: 1 genome(s) have been classified using the ANI pre-screening step.
    [2023-02-22 16:20:11] INFO: Done.
    [2023-02-22 16:20:11] INFO: 1 genome(s) have been classified using the ANI pre-screening step.
    [2023-02-22 16:20:11] INFO: Done.
    [2023-02-22 16:20:11] INFO: Identifying markers in 2 genomes with 20 threads.
    [2023-02-22 16:20:11] TASK: Running Prodigal V2.6.3 to identify genes.
    [2023-02-22 16:20:12] INFO: Completed 2 genomes in 0.22 seconds (9.07 genomes/second).
    [2023-02-22 16:20:12] WARNING: Prodigal skipped 2 genomes due to pre-existing data, see warnings.log
    [2023-02-22 16:20:12] TASK: Identifying TIGRFAM protein families.
    [2023-02-22 16:20:12] INFO: Completed 2 genomes in 0.03 seconds (65.39 genomes/second).
    [2023-02-22 16:20:12] WARNING: TIGRFAM skipped 2 genomes due to pre-existing data, see warnings.log
    [2023-02-22 16:20:12] TASK: Identifying Pfam protein families.
    [2023-02-22 16:20:12] INFO: Completed 2 genomes in 0.03 seconds (68.36 genomes/second).
    [2023-02-22 16:20:12] WARNING: Pfam skipped 2 genomes due to pre-existing data, see warnings.log
    [2023-02-22 16:20:12] INFO: Annotations done using HMMER 3.1b2 (February 2015).
    [2023-02-22 16:20:12] TASK: Summarising identified marker genes.
    [2023-02-22 16:20:12] INFO: Completed 2 genomes in 0.06 seconds (32.55 genomes/second).
    [2023-02-22 16:20:12] INFO: Done.
    [2023-02-22 16:20:12] INFO: Aligning markers in 2 genomes with 20 CPUs.
    [2023-02-22 16:20:12] INFO: Processing 2 genomes identified as bacterial.
    [2023-02-22 16:20:21] INFO: Read concatenated alignment for 62,291 GTDB genomes.
    [2023-02-22 16:20:21] TASK: Generating concatenated alignment for each marker.
    [2023-02-22 16:20:22] INFO: Completed 2 genomes in 0.03 seconds (79.85 genomes/second).
    [2023-02-22 16:20:23] TASK: Aligning 100 identified markers using hmmalign 3.1b2 (February 2015).
    [2023-02-22 16:20:25] INFO: Completed 100 markers in 1.06 seconds (93.94 markers/second).
    [2023-02-22 16:20:25] TASK: Masking columns of bacterial multiple sequence alignment using canonical mask.
    [2023-02-22 16:22:21] INFO: Completed 62,293 sequences in 1.93 minutes (32,233.24 sequences/minute).
    [2023-02-22 16:22:21] INFO: Masked bacterial alignment from 41,084 to 5,036 AAs.
    [2023-02-22 16:22:21] INFO: 0 bacterial user genomes have amino acids in <10.0% of columns in filtered MSA.
    [2023-02-22 16:22:22] INFO: Creating concatenated alignment for 62,293 bacterial GTDB and user genomes.
    [2023-02-22 16:22:46] INFO: Creating concatenated alignment for 2 bacterial user genomes.
    [2023-02-22 16:22:46] INFO: Done.
    [2023-02-22 16:22:47] TASK: Placing 2 bacterial genomes into backbone reference tree with pplacer using 20 CPUs (be patient).
    [2023-02-22 16:22:47] INFO: pplacer version: v1.1.alpha19-0-g807f6f3
    [2023-02-22 16:25:01] INFO: Calculating RED values based on reference tree.
    [2023-02-22 16:25:02] INFO: 2 out of 2 have an class assignments. Those genomes will be reclassified.
    [2023-02-22 16:25:02] TASK: Placing 1 bacterial genomes into class-level reference tree 6 (1/2) with pplacer using 20 CPUs (be patient).
    [2023-02-22 16:29:46] INFO: Calculating RED values based on reference tree.
    [2023-02-22 16:29:48] TASK: Traversing tree to determine classification method.
    [2023-02-22 16:29:48] INFO: Completed 1 genome in 0.00 seconds (2,391.28 genomes/second).
    [2023-02-22 16:29:48] TASK: Calculating average nucleotide identity using FastANI (v1.3).
    [2023-02-22 16:29:50] INFO: Completed 34 comparisons in 1.53 seconds (22.22 comparisons/second).
    [2023-02-22 16:29:50] INFO: 0 genome(s) have been classified using FastANI and pplacer.
    [2023-02-22 16:29:50] TASK: Placing 1 bacterial genomes into class-level reference tree 5 (2/2) with pplacer using 20 CPUs (be patient).
    [2023-02-22 16:33:17] INFO: Calculating RED values based on reference tree.
    [2023-02-22 16:33:19] TASK: Traversing tree to determine classification method.
    [2023-02-22 16:33:19] INFO: Completed 1 genome in 0.06 seconds (17.02 genomes/second).
    [2023-02-22 16:33:19] INFO: 0 genome(s) have been classified using FastANI and pplacer.
    [2023-02-22 16:33:19] WARNING: 1 of 3 genome has a warning (see summary file).
    [2023-02-22 16:33:19] INFO: 0 genome(s) have been classified using FastANI and pplacer.
    [2023-02-22 16:33:19] WARNING: 1 of 3 genome has a warning (see summary file).
    [2023-02-22 16:33:19] INFO: Note that Tk classification mode is insufficient for publication of new taxonomic designations. New designations should be based on one or more de novo trees, an example of which can be produced by Tk in de novo mode.
    [2023-02-22 16:33:19] INFO: Done.
