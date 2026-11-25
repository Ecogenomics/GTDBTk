.. _commands/classify:

classify
========

Determine taxonomic classification of genomes.

Arguments
---------

.. argparse::
   :module: gtdbtk.cli
   :func: get_main_parser
   :prog: gtdbtk
   :path: classify
   :nodefaultconst:


Files output
------------

* classify
    * :ref:`[prefix].[domain].summary.tsv <files/summary.tsv>`
    * :ref:`[prefix].[domain].classify.tree <files/classify.tree>`
    * intermediate_results
        * :ref:`[prefix].[domain].classification_gplace.tsv <files/classification_gplace.tsv>`
        * :ref:`[prefix].[domain].prescreened.msa.fasta <files/msa.fasta>`
        * :ref:`[prefix].[domain].red_dictionary.tsv <files/red_dictionary.tsv>`
        * gplace
            * :ref:`gplace.[domain].json <files/gplace.domain.json>`
            * :ref:`gplace.[domain].out <files/gplace.domain.out>`
* ani_screen
* :ref:`[prefix].[domain].summary.tsv <files/summary.tsv>`
* :ref:`[prefix].log <files/gtdbtk.log>`
* :ref:`[prefix].json <files/gtdbtk.json>`
* :ref:`[prefix].warnings.log <files/gtdbtk.warnings.log>`


Example
-------


Input
^^^^^


.. code-block:: bash

     gtdbtk classify --batchfile genomes/3_batchfile.tsv --align_dir 3_align/ --out_dir 3_classify --cpus 50



Output
^^^^^^


.. code-block:: text

    [2025-08-05 17:41:11] INFO: GTDB-Tk v2.5.0
    [2025-08-05 17:41:11] INFO: gtdbtk classify --batchfile genomes/3_batchfile.tsv --align_dir 3_align/ --out_dir 3_classify --cpus 50
    [2025-08-05 17:41:11] INFO: Using GTDB-Tk reference data version r226: /srv/db/gtdbtk/official/release226
    [2025-08-05 17:41:13] INFO: Loading reference genomes.
    [2025-08-05 17:41:13] INFO: Calculating all vs all ANI with skani v0.2.1.
    [2025-08-05 17:41:14] INFO: Sketching genomes
    [2025-08-05 17:42:21] INFO: Sketches done: 1min 6secs
    [2025-08-05 17:42:22] INFO: Running comparisons
    [2025-08-05 17:42:51] INFO: Comparisons finished, capturing results.
    [2025-08-05 17:44:45] INFO: 0 genome(s) have been classified using the ANI pre-screening step.
    [2025-08-05 17:44:45] TASK: Placing 3 bacterial genomes into reference tree with gplace using 50 CPUs (be patient).
    [2025-08-05 17:44:45] INFO: gplace version: 0.1.1
    [2025-08-05 17:48:20] INFO: Calculating RED values based on reference tree.
    [2025-08-05 17:48:26] TASK: Traversing tree to determine classification method.
    [2025-08-05 17:48:26] INFO: Completed 3 genomes in 0.00 seconds (6,043.67 genomes/second).
    [2025-08-05 17:48:26] TASK: Calculating average nucleotide identity using skani (v0.2.1).
    [2025-08-05 17:48:27] INFO: Completed 13 comparisons in 1.18 seconds (11.02 comparisons/second).
    [2025-08-05 17:48:27] INFO: 0 genome(s) have been classified using skani and gplace.
    [2025-08-05 17:48:27] INFO: Note that Tk classification mode is insufficient for publication of new taxonomic designations. New designations should be based on one or more de novo trees, an example of which can be produced by Tk in de novo mode.
    [2025-08-05 17:48:27] INFO: Done.