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
    * intermediate_results
        * mash
            * :ref:`[prefix].mash_distances.tsv <files/mash_distances.msh>`
            * :ref:`[prefix].user_query_sketch.msh <files/user_query_sketch.msh>`
* :ref:`[prefix].[domain].summary.tsv <files/summary.tsv>`
* :ref:`[prefix].log <files/gtdbtk.log>`
* :ref:`[prefix].json <files/gtdbtk.json>`
* :ref:`[prefix].warnings.log <files/gtdbtk.warnings.log>`


Example
-------


Input
^^^^^


.. code-block:: bash

     gtdbtk classify --align_dir align_3lines/ --batchfile 3lines_batchfile.tsv --out_dir 3classify_ani --mash_db mash_db_dir/ --cpus 20



Output
^^^^^^


.. code-block:: text

    [2023-02-15 08:37:11] INFO: GTDB-Tk v2.2.2
    [2023-02-15 08:37:11] INFO: gtdbtk classify --align_dir align_3lines/ --batchfile 3lines_batchfile.tsv --out_dir 3classify_ani --mash_db mash_db_dir/ --cpus 20
    [2023-02-15 08:37:11] INFO: Using GTDB-Tk reference data version r207: /srv/projects/gtdbtk/test_new_features/release207_v2/
    [2023-02-15 08:37:12] INFO: Loading reference genomes.
    [2023-02-15 08:37:13] INFO: Using Mash version 2.2.2
    [2023-02-15 08:37:13] INFO: Loading data from existing Mash sketch file: 3classify_ani/classify/ani_screen/intermediate_results/mash/gtdbtk.user_query_sketch.msh
    [2023-02-15 08:37:13] INFO: Loading data from existing Mash sketch file: mash_db_dir/gtdb_ref_sketch.msh
    [2023-02-15 08:37:16] INFO: Calculating Mash distances.
    [2023-02-15 08:37:20] INFO: Calculating ANI with FastANI v1.3.
    [2023-02-15 08:37:21] INFO: Completed 12 comparisons in 0.62 seconds (19.21 comparisons/second).
    [2023-02-15 08:37:21] INFO: 1 genome(s) have been classified using the ANI pre-screening step.
    [2023-02-15 08:37:21] TASK: Placing 2 bacterial genomes into backbone reference tree with pplacer using 20 CPUs (be patient).
    [2023-02-15 08:37:21] INFO: pplacer version: v1.1.alpha19-0-g807f6f3
    [2023-02-15 08:39:24] INFO: Calculating RED values based on reference tree.
    [2023-02-15 08:39:25] INFO: 2 out of 2 have an class assignments. Those genomes will be reclassified.
    [2023-02-15 08:39:25] TASK: Placing 1 bacterial genomes into class-level reference tree 6 (1/2) with pplacer using 20 CPUs (be patient).
    [2023-02-15 08:43:39] INFO: Calculating RED values based on reference tree.
    [2023-02-15 08:43:42] TASK: Traversing tree to determine classification method.
    [2023-02-15 08:43:42] INFO: Completed 1 genome in 0.00 seconds (2,451.38 genomes/second).
    [2023-02-15 08:43:42] TASK: Calculating average nucleotide identity using FastANI (v1.3).
    [2023-02-15 08:43:43] INFO: Completed 34 comparisons in 0.90 seconds (37.77 comparisons/second).
    [2023-02-15 08:43:43] INFO: 0 genome(s) have been classified using FastANI and pplacer.
    [2023-02-15 08:43:43] TASK: Placing 1 bacterial genomes into class-level reference tree 5 (2/2) with pplacer using 20 CPUs (be patient).
    [2023-02-15 08:46:38] INFO: Calculating RED values based on reference tree.
    [2023-02-15 08:46:40] TASK: Traversing tree to determine classification method.
    [2023-02-15 08:46:40] INFO: Completed 1 genome in 0.05 seconds (20.80 genomes/second).
    [2023-02-15 08:46:40] INFO: 0 genome(s) have been classified using FastANI and pplacer.
    [2023-02-15 08:46:41] WARNING: 1 of 3 genome has a warning (see summary file).
    [2023-02-15 08:46:41] INFO: Note that Tk classification mode is insufficient for publication of new taxonomic designations. New designations should be based on one or more de novo trees, an example of which can be produced by Tk in de novo mode.
    [2023-02-15 08:46:41] INFO: Done.