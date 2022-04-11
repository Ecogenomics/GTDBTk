.. _commands/identify:

identify
========

Identify marker genes in genome(s). The following heuristic is used to establish the translation table used by Prodigal: use table 11 unless the coding density using table 4 is 5% higher than when using table 11 and the coding density under table 4 is >70%. Distinguishing between tables 4 and 25 is challenging so GTDB-Tk does not attempt to distinguish between these two tables. If you know the correct translation table for your genomes this can be provided to GTDB-Tk in the `--batchfile`.

Arguments
---------

.. argparse::
   :module: gtdbtk.cli
   :func: get_main_parser
   :prog: gtdbtk
   :path: identify
   :nodefaultconst:


## Files output

* :ref:`[prefix].log <files/gtdbtk.log>`
* :ref:`[prefix].warnings.log <files/gtdbtk.warnings.log>`
* identify/
    * :ref:`[prefix].[domain].markers_summary.tsv <files/markers_summary.tsv>`
    * :ref:`[prefix].translation_table_summary.tsv <files/translation_table_summary.tsv>`
* identify/intermediate_results/marker_genes/[genome_id]/
    * :ref:`[genome_id]_pfam_tophit.tsv <files/pfam_tophit.tsv>`
    * :ref:`[genome_id]_pfam.tsv <files/pfam.tsv>`
    * :ref:`[genome_id]_protein.faa <files/protein.faa>`
    * :ref:`[genome_id]_protein.fna <files/protein.fna>`
    * :ref:`[genome_id]_protein.gff <files/protein.gff>`
    * :ref:`[genome_id]_tigrfam.out <files/tigrfam.out>`
    * :ref:`[genome_id]_tigrfam_tophit.tsv <files/tigrfam_tophit.tsv>`
    * :ref:`[genome_id]_tigrfam.tsv <files/tigrfam.tsv>`
    * :ref:`prodigal_translation_table.tsv <files/prodigal_translation_table.tsv>`

Example
-------


Input
^^^^^


.. code-block:: bash

    gtdbtk identify --genome_dir genomes/ --out_dir identify_output --cpus 3



Output
^^^^^^

.. code-block:: text

    [2022-04-11 11:48:59] INFO: GTDB-Tk v2.0.0
    [2022-04-11 11:48:59] INFO: gtdbtk identify --genome_dir /tmp/gtdbtk/genomes --out_dir /tmp/gtdbtk/identify --extension gz --cpus 2
    [2022-04-11 11:48:59] INFO: Using GTDB-Tk reference data version r207: /srv/db/gtdbtk/official/release207
    [2022-04-11 11:48:59] INFO: Identifying markers in 2 genomes with 2 threads.
    [2022-04-11 11:48:59] TASK: Running Prodigal V2.6.3 to identify genes.
    [2022-04-11 11:49:10] INFO: Completed 2 genomes in 10.94 seconds (5.47 seconds/genome).
    [2022-04-11 11:49:10] TASK: Identifying TIGRFAM protein families.
    [2022-04-11 11:49:16] INFO: Completed 2 genomes in 5.78 seconds (2.89 seconds/genome).
    [2022-04-11 11:49:16] TASK: Identifying Pfam protein families.
    [2022-04-11 11:49:16] INFO: Completed 2 genomes in 0.42 seconds (4.81 genomes/second).
    [2022-04-11 11:49:16] INFO: Annotations done using HMMER 3.1b2 (February 2015).
    [2022-04-11 11:49:16] TASK: Summarising identified marker genes.
    [2022-04-11 11:49:16] INFO: Completed 2 genomes in 0.05 seconds (40.91 genomes/second).
    [2022-04-11 11:49:16] INFO: Done.

