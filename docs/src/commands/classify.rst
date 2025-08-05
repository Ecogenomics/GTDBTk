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
    * :ref:`[prefix].backbone.[domain].classify.tree <files/classify.tree>`
    * :ref:`[prefix].[domain].tree.mapping.tsv <files/tree.mapping.tsv>`
    * :ref:`[prefix].[domain].classify.tree.[index].tree <files/classify.tree>`
    * intermediate_results
        * :ref:`[prefix].[domain].backbone.classification_pplacer.tsv <files/classification_pplacer.tsv>`
        * :ref:`[prefix].[domain].class_level.classification_pplacer_tree_[index].tsv <files/classification_pplacer.tsv>`
        * :ref:`[prefix].[domain].prescreened.msa.fasta <files/msa.fasta>`
        * :ref:`[prefix].[domain].red_dictionary.tsv <files/red_dictionary.tsv>`
        * pplacer
            * :ref:`pplacer.backbone.[domain].json <files/pplacer.domain.json>`
            * :ref:`pplacer.backbone.[domain].out <files/pplacer.domain.out>`
            * tree_[index]
                * :ref:`[prefix].[domain].user_msa.fasta <files/user_msa.fasta>`
                * :ref:`pplacer.class_level.[domain].out <files/pplacer.domain.out>`
                * :ref:`pplacer.class_level.[domain].json <files/pplacer.domain.json>`
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
    [2025-08-05 17:44:45] TASK: Placing 3 bacterial genomes into backbone reference tree with pplacer using 50 CPUs (be patient).
    [2025-08-05 17:44:45] INFO: pplacer version: v1.1.alpha19-0-g807f6f3
    [2025-08-05 17:48:20] INFO: Calculating RED values based on reference tree.
    [2025-08-05 17:48:21] INFO: 3 out of 3 have an class assignments. Those genomes will be reclassified.
    [2025-08-05 17:48:21] TASK: Placing 2 bacterial genomes into class-level reference tree 7 (1/2) with pplacer using 50 CPUs (be patient).
    [2025-08-05 18:00:49] INFO: Calculating RED values based on reference tree.
    [2025-08-05 18:00:54] TASK: Traversing tree to determine classification method.
    [2025-08-05 18:00:54] INFO: Completed 2 genomes in 0.00 seconds (6,043.67 genomes/second).
    [2025-08-05 18:00:54] TASK: Calculating average nucleotide identity using skani (v0.2.1).
    [2025-08-05 18:00:55] INFO: Completed 13 comparisons in 1.18 seconds (11.02 comparisons/second).
    [2025-08-05 18:00:55] INFO: 0 genome(s) have been classified using skani and pplacer.
    [2025-08-05 18:00:56] TASK: Placing 1 bacterial genomes into class-level reference tree 2 (2/2) with pplacer using 50 CPUs (be patient).
    [2025-08-05 18:12:59] INFO: Calculating RED values based on reference tree.
    [2025-08-05 18:13:05] TASK: Traversing tree to determine classification method.
    [2025-08-05 18:13:05] INFO: Completed 1 genome in 0.00 seconds (1,681.08 genomes/second).
    [2025-08-05 18:13:05] INFO: 0 genome(s) have been classified using skani and pplacer.
    [2025-08-05 18:13:05] INFO: Note that Tk classification mode is insufficient for publication of new taxonomic designations. New designations should be based on one or more de novo trees, an example of which can be produced by Tk in de novo mode.
    [2025-08-05 18:13:05] INFO: Done.