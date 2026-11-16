.. _files/gplace.domain.out:

gplace.<domain>.out
===================

The log written by gPlace while placing the genomes (its warnings are also copied to the GTDB-Tk log).

Produced by
-----------
 * :ref:`commands/classify`
 * :ref:`commands/classify_wf`


Example
-------

.. code-block:: text

    [2026-11-13 10:00:00] INFO: gPlace v0.1.1
    [2026-11-13 10:00:00] INFO: gplace place -m wag -j 3 -c /release232/pplacer/gtdb_r232_ar53.refpkg -o classify_output/classify/intermediate_results/gplace/gplace.ar53.json align_output/align/gtdbtk.ar53.user_msa.fasta.gz
    [2026-11-13 10:00:00] INFO: --model: no effect, the reference package's phylo_model sets the model (as in pplacer).
    [2026-11-13 10:00:00] INFO: Read 3 sequences from align_output/align/gtdbtk.ar53.user_msa.fasta.gz
    [2026-11-13 10:00:01] INFO: Model: WAG + CAT (20 rates), model frequencies (protein)
    [2026-11-13 10:00:01] INFO: Placing 3 queries on /release232/pplacer/gtdb_r232_ar53.refpkg (…  nodes, … sites, 3 CPUs; estimated peak memory … GiB).
    [2026-11-13 10:00:20] INFO: Cached likelihood vectors for … internal nodes (… GB).
    [2026-11-13 10:00:25] INFO: Scored 3 queries against … edges.
    [2026-11-13 10:00:31] INFO: Placed 3/3 queries.
    [2026-11-13 10:00:31] INFO: Wrote 3 placements (3 records) to classify_output/classify/intermediate_results/gplace/gplace.ar53.json (31.0 s).
