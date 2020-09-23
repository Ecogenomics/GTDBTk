.. _commands/classify_wf:

classify_wf
===========

Classify workflow
-----------------

For arguments and output files, see each of the individual steps:

* :ref:`commands/infer`
* :ref:`commands/align`
* :ref:`commands/classify`

The classify workflow consists of three steps: ``identify``, ``align``, and ``classify``.

The ``identify`` step calls genes using `Prodigal <http://compbio.ornl.gov/prodigal/>`_,
and uses HMM models and the `HMMER <http://hmmer.org/>`_ package to identify the
120 bacterial and 122 archaeal marker genes used for phylogenetic inference
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
   :module: gtdbtk.argparse
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

    [2020-04-14 09:46:06] INFO: GTDB-Tk v1.1.0
    [2020-04-14 09:46:06] INFO: gtdbtk classify_wf --genome_dir genomes/ --out_dir classify_wf_out --cpus 3
    [2020-04-14 09:46:06] INFO: Using GTDB-Tk reference data version r89: /release89
    [2020-04-14 09:46:06] INFO: Identifying markers in 3 genomes with 3 threads.
    [2020-04-14 09:46:06] INFO: Running Prodigal V2.6.3 to identify genes.
    ==> Finished processing 3 of 3 (100.0%) genomes.
    [2020-04-14 09:46:24] INFO: Identifying TIGRFAM protein families.
    ==> Finished processing 3 of 3 (100.0%) genomes.
    [2020-04-14 09:46:33] INFO: Identifying Pfam protein families.
    ==> Finished processing 3 of 3 (100.0%) genomes.
    [2020-04-14 09:46:35] INFO: Annotations done using HMMER 3.1b2 (February 2015)
    [2020-04-14 09:46:35] INFO: Done.
    [2020-04-14 09:46:35] INFO: Aligning markers in 3 genomes with 3 threads.
    [2020-04-14 09:46:35] INFO: Processing 3 genomes identified as archaeal.
    [2020-04-14 09:46:35] INFO: Read concatenated alignment for 1248 GTDB genomes.
    ==> Finished aligning 3 of 3 (100.0%) genomes.
    [2020-04-14 09:46:40] INFO: Masking columns of multiple sequence alignment using canonical mask.
    [2020-04-14 09:46:42] INFO: Masked alignment from 32675 to 5124 AAs.
    [2020-04-14 09:46:42] INFO: 0 user genomes have amino acids in <10.0% of columns in filtered MSA.
    [2020-04-14 09:46:42] INFO: Creating concatenated alignment for 1251 GTDB and user genomes.
    [2020-04-14 09:46:42] INFO: Creating concatenated alignment for 3 user genomes.
    [2020-04-14 09:46:42] INFO: Done.
    [2020-04-14 09:46:42] INFO: Placing 3 archaeal genomes into reference tree with pplacer using 3 cpus (be patient).
    Placing genomes |##################################################| 3/3 (100.00%)
    [2020-04-14 09:47:41] INFO: pplacer version: v1.1.alpha19-0-g807f6f3
    [2020-04-14 09:47:41] INFO: Calculating average nucleotide identity using FastANI.
    [2020-04-14 09:47:41] INFO: fastANI version: 1.3
    ==> Processing 24 of 24 (100.0%) comparisons.
    [2020-04-14 09:47:46] INFO: 3 genome(s) have been classified using FastANI and pplacer.
    [2020-04-14 09:47:46] INFO: Calculating RED values based on reference tree.
    [2020-04-14 09:47:47] INFO: Done.

