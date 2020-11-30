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

* :ref:`[prefix].[domain].markers_summary.tsv <files/markers_summary.tsv>`
* :ref:`[prefix].log <files/gtdbtk.log>`
* :ref:`[prefix].translation_table_summary.tsv <files/translation_table_summary.tsv>`
* :ref:`[prefix].warnings.log <files/gtdbtk.warnings.log>`
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
    
    [2020-04-14 08:51:00] INFO: GTDB-Tk v1.1.0
    [2020-04-14 08:51:00] INFO: gtdbtk identify --genome_dir genomes/ --out_dir identify_output --cpus 3
    [2020-04-14 08:51:00] INFO: Using GTDB-Tk reference data version r89: /release89
    [2020-04-14 08:51:00] INFO: Identifying markers in 3 genomes with 3 threads.
    [2020-04-14 08:51:00] INFO: Running Prodigal V2.6.3 to identify genes.
    ==> Finished processing 3 of 3 (100.0%) genomes.
    [2020-04-14 08:51:18] INFO: Identifying TIGRFAM protein families.
    ==> Finished processing 3 of 3 (100.0%) genomes.
    [2020-04-14 08:51:27] INFO: Identifying Pfam protein families.
    ==> Finished processing 3 of 3 (100.0%) genomes.
    [2020-04-14 08:51:29] INFO: Annotations done using HMMER 3.1b2 (February 2015)
    [2020-04-14 08:51:29] INFO: Done.
    
