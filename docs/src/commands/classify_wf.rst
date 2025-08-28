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

The ``ani_screen`` step compares user genomes against a `skani <https://www.nature.com/articles/s41592-023-02018-3>`_ database composed of all GTDB representative genomes.
User genomes classified with skani are not run through the rest of the pipeline (``identify``, ``align``, ``classify``)
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

    gtdbtk classify_wf --batchfile genomes/3_batchfile.tsv --out_dir classify_wf_3_genomes --cpus 20



Output
^^^^^^


.. code-block:: text

    [2025-08-05 20:31:02] INFO: GTDB-Tk v2.5.0
    [2025-08-05 20:31:02] INFO: gtdbtk classify_wf --batchfile genomes/3_batchfile.tsv --out_dir classify_wf_3_genomes --cpus 20
    [2025-08-05 20:31:02] INFO: Using GTDB-Tk reference data version r226: /srv/db/gtdbtk/official/release226
    [2025-08-05 20:31:02] INFO: Loading reference genomes.
    [2025-08-05 20:31:03] INFO: Calculating all vs all ANI with skani v0.2.1.
    [2025-08-05 20:31:04] INFO: Sketching genomes
    [2025-08-05 20:33:38] INFO: Sketches done: 2min 34secs
    [2025-08-05 20:33:38] INFO: Running comparisons
    [2025-08-05 20:33:57] INFO: Comparisons finished, capturing results.
    [2025-08-05 20:35:27] INFO: 0 genome(s) have been classified using the ANI pre-screening step.
    [2025-08-05 20:35:27] INFO: Done.
    [2025-08-05 20:35:27] INFO: Identifying markers in 3 genomes with 20 threads.
    [2025-08-05 20:35:27] TASK: Running Prodigal V2.6.3 to identify genes.
    [2025-08-05 20:36:09] INFO: Completed 3 genomes in 41.78 seconds (13.93 seconds/genome).
    [2025-08-05 20:36:09] TASK: Identifying TIGRFAM protein families.
    [2025-08-05 20:36:14] INFO: Completed 3 genomes in 5.57 seconds (1.86 seconds/genome).
    [2025-08-05 20:36:14] TASK: Identifying Pfam protein families.
    [2025-08-05 20:36:15] INFO: Completed 3 genomes in 0.70 seconds (4.28 genomes/second).
    [2025-08-05 20:36:15] INFO: Annotations done using HMMER 3.1b2 (February 2015).
    [2025-08-05 20:36:15] TASK: Summarising identified marker genes.
    [2025-08-05 20:36:15] INFO: Completed 3 genomes in 0.17 seconds (17.58 genomes/second).
    [2025-08-05 20:36:15] INFO: Done.
    [2025-08-05 20:36:15] INFO: Aligning markers in 3 genomes with 20 CPUs.
    [2025-08-05 20:36:16] INFO: Processing 3 genomes identified as bacterial.
    [2025-08-05 20:36:27] INFO: Read concatenated alignment for 136,646 GTDB genomes.
    [2025-08-05 20:36:27] TASK: Generating concatenated alignment for each marker.
    [2025-08-05 20:36:33] INFO: Completed 3 genomes in 0.07 seconds (45.19 genomes/second).
    [2025-08-05 20:36:34] TASK: Aligning 119 identified markers using hmmalign 3.1b2 (February 2015).
    [2025-08-05 20:36:46] INFO: Completed 119 markers in 5.66 seconds (21.02 markers/second).
    [2025-08-05 20:36:46] TASK: Masking columns of bacterial multiple sequence alignment using canonical mask.
    [2025-08-05 20:40:54] INFO: Completed 136,649 sequences in 4.12 minutes (33,135.45 sequences/minute).
    [2025-08-05 20:40:54] INFO: Masked bacterial alignment from 41,084 to 5,036 AAs.
    [2025-08-05 20:40:54] INFO: 0 bacterial user genomes have amino acids in <10.0% of columns in filtered MSA.
    [2025-08-05 20:40:54] INFO: Creating concatenated alignment for 136,649 bacterial GTDB and user genomes.
    [2025-08-05 20:41:36] INFO: Creating concatenated alignment for 3 bacterial user genomes.
    [2025-08-05 20:41:38] INFO: Done.
    [2025-08-05 20:41:39] TASK: Placing 3 bacterial genomes into backbone reference tree with pplacer using 20 CPUs (be patient).
    [2025-08-05 20:41:39] INFO: pplacer version: v1.1.alpha19-0-g807f6f3
    [2025-08-05 20:44:58] INFO: Calculating RED values based on reference tree.
    [2025-08-05 20:44:59] INFO: 3 out of 3 have an class assignments. Those genomes will be reclassified.
    [2025-08-05 20:44:59] TASK: Placing 2 bacterial genomes into class-level reference tree 7 (1/2) with pplacer using 20 CPUs (be patient).
    [2025-08-05 20:56:48] INFO: Calculating RED values based on reference tree.
    [2025-08-05 20:56:52] TASK: Traversing tree to determine classification method.
    [2025-08-05 20:56:52] INFO: Completed 2 genomes in 0.00 seconds (6,026.30 genomes/second).
    [2025-08-05 20:56:52] TASK: Calculating average nucleotide identity using skani (v0.2.1).
    [2025-08-05 20:56:53] INFO: Completed 13 comparisons in 0.53 seconds (24.35 comparisons/second).
    [2025-08-05 20:56:53] INFO: 0 genome(s) have been classified using skani and pplacer.
    [2025-08-05 20:56:53] TASK: Placing 1 bacterial genomes into class-level reference tree 2 (2/2) with pplacer using 20 CPUs (be patient).
    [2025-08-05 21:08:19] INFO: Calculating RED values based on reference tree.
    [2025-08-05 21:08:25] TASK: Traversing tree to determine classification method.
    [2025-08-05 21:08:25] INFO: Completed 1 genome in 0.00 seconds (1,531.89 genomes/second).
    [2025-08-05 21:08:25] INFO: 0 genome(s) have been classified using skani and pplacer.
    [2025-08-05 21:08:25] INFO: Note that Tk classification mode is insufficient for publication of new taxonomic designations. New designations should be based on one or more de novo trees, an example of which can be produced by Tk in de novo mode.
    [2025-08-05 21:08:25] INFO: Done.
    [2025-08-05 21:08:25] INFO: Removing intermediate files.
    [2025-08-05 21:08:25] INFO: Intermediate files removed.
    [2025-08-05 21:08:25] INFO: Done.
