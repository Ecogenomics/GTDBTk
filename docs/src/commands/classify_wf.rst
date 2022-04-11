.. _commands/classify_wf:

classify_wf
===========

Classify workflow
-----------------

For arguments and output files, see each of the individual steps:

* :ref:`commands/identify`
* :ref:`commands/align`
* :ref:`commands/classify`

The classify workflow consists of three steps: ``identify``, ``align``, and ``classify``.

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

    [2022-04-11 12:48:53] INFO: GTDB-Tk v2.0.0
    [2022-04-11 12:48:53] INFO: gtdbtk classify_wf --genome_dir genomes/ --out_dir classify_wf_out --cpus 3 -x gz
    [2022-04-11 12:48:53] INFO: Using GTDB-Tk reference data version r207: /srv/db/gtdbtk/official/release207
    [2022-04-11 12:48:53] INFO: Identifying markers in 3 genomes with 3 threads.
    [2022-04-11 12:48:53] TASK: Running Prodigal V2.6.3 to identify genes.
    [2022-04-11 12:49:04] INFO: Completed 3 genomes in 10.96 seconds (3.65 seconds/genome).
    [2022-04-11 12:49:04] TASK: Identifying TIGRFAM protein families.
    [2022-04-11 12:49:10] INFO: Completed 3 genomes in 5.88 seconds (1.96 seconds/genome).
    [2022-04-11 12:49:10] TASK: Identifying Pfam protein families.
    [2022-04-11 12:49:10] INFO: Completed 3 genomes in 0.41 seconds (7.30 genomes/second).
    [2022-04-11 12:49:10] INFO: Annotations done using HMMER 3.1b2 (February 2015).
    [2022-04-11 12:49:10] TASK: Summarising identified marker genes.
    [2022-04-11 12:49:11] INFO: Completed 3 genomes in 0.07 seconds (40.18 genomes/second).
    [2022-04-11 12:49:11] INFO: Done.
    [2022-04-11 12:49:11] INFO: Aligning markers in 3 genomes with 3 CPUs.
    [2022-04-11 12:49:11] INFO: Processing 3 genomes identified as archaeal.
    [2022-04-11 12:49:11] INFO: Read concatenated alignment for 3,412 GTDB genomes.
    [2022-04-11 12:49:11] TASK: Generating concatenated alignment for each marker.
    [2022-04-11 12:49:11] INFO: Completed 3 genomes in 0.02 seconds (167.25 genomes/second).
    [2022-04-11 12:49:11] TASK: Aligning 52 identified markers using hmmalign 3.1b2 (February 2015).
    [2022-04-11 12:49:11] INFO: Completed 52 markers in 0.54 seconds (96.16 markers/second).
    [2022-04-11 12:49:11] TASK: Masking columns of archaeal multiple sequence alignment using canonical mask.
    [2022-04-11 12:49:16] INFO: Completed 3,415 sequences in 4.15 seconds (822.38 sequences/second).
    [2022-04-11 12:49:16] INFO: Masked archaeal alignment from 13,540 to 10,153 AAs.
    [2022-04-11 12:49:16] INFO: 0 archaeal user genomes have amino acids in <10.0% of columns in filtered MSA.
    [2022-04-11 12:49:16] INFO: Creating concatenated alignment for 3,415 archaeal GTDB and user genomes.
    [2022-04-11 12:49:18] INFO: Creating concatenated alignment for 3 archaeal user genomes.
    [2022-04-11 12:49:18] INFO: Done.
    [2022-04-11 12:49:18] TASK: Placing 3 archaeal genomes into reference tree with pplacer using 3 CPUs (be patient).
    [2022-04-11 12:49:18] INFO: pplacer version: v1.1.alpha19-0-g807f6f3
    [2022-04-11 12:54:22] INFO: Calculating RED values based on reference tree.
    [2022-04-11 12:54:23] TASK: Traversing tree to determine classification method.
    [2022-04-11 12:54:23] INFO: Completed 3 genomes in 0.00 seconds (23,563.51 genomes/second).
    [2022-04-11 12:54:23] TASK: Calculating average nucleotide identity using FastANI (v1.32).
    [2022-04-11 12:54:25] INFO: Completed 6 comparisons in 1.96 seconds (3.06 comparisons/second).
    [2022-04-11 12:54:25] INFO: 3 genome(s) have been classified using FastANI and pplacer.
    [2022-04-11 12:54:25] INFO: Note that Tk classification mode is insufficient for publication of new taxonomic designations. New designations should be based on one or more de novo trees, an example of which can be produced by Tk in de novo mode.
    [2022-04-11 12:54:25] INFO: Done.
    [2022-04-11 12:54:25] INFO: Removing intermediate files.
    [2022-04-11 12:54:25] INFO: Intermediate files removed.
    [2022-04-11 12:54:25] INFO: Done.
