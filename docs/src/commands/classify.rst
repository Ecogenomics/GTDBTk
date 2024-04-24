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

    [2024-03-27 15:46:26] INFO: GTDB-Tk v2.3.2
    [2024-03-27 15:46:26] INFO: gtdbtk classify --align_dir 500_align --out_dir 500_classify --mash_db mash_db.msh --cpus 90 --batchfile genomes/500_batchfile.tsv
    [2024-03-27 15:46:26] INFO: Using GTDB-Tk reference data version r214: /srv/db/gtdbtk/official/release214_skani/release214
    [2024-03-27 15:46:27] WARNING: Setting pplacer CPUs to 64, as pplacer is known to hang if >64 are used. You can override this using: --pplacer_cpus
    [2024-03-27 15:46:27] INFO: Loading reference genomes.
    [2024-03-27 15:46:27] INFO: Using Mash version 2.2.2
    [2024-03-27 15:46:27] INFO: Creating Mash sketch file: 500_classify/classify/ani_screen/intermediate_results/mash/gtdbtk.user_query_sketch.msh
    [2024-03-27 15:46:29] INFO: Completed 500 genomes in 1.74 seconds (286.68 genomes/second).
    [2024-03-27 15:46:29] INFO: Loading data from existing Mash sketch file: mash_db.msh
    [2024-03-27 15:46:32] INFO: Calculating Mash distances.
    [2024-03-27 15:47:17] INFO: Calculating ANI with skani v0.2.1.
    [2024-03-27 15:47:28] INFO: Completed 4,383 comparisons in 10.51 seconds (417.14 comparisons/second).
    [2024-03-27 15:47:30] INFO: 357 genome(s) have been classified using the ANI pre-screening step.
    [2024-03-27 15:47:30] TASK: Placing 143 bacterial genomes into backbone reference tree with pplacer using 64 CPUs (be patient).
    [2024-03-27 15:47:30] INFO: pplacer version: v1.1.alpha19-0-g807f6f3
    [2024-03-27 15:50:19] INFO: Calculating RED values based on reference tree.
    [2024-03-27 15:50:20] INFO: 143 out of 143 have an class assignments. Those genomes will be reclassified.
    [2024-03-27 15:50:20] TASK: Placing 30 bacterial genomes into class-level reference tree 3 (1/8) with pplacer using 64 CPUs (be patient).
    [2024-03-27 15:58:09] INFO: Calculating RED values based on reference tree.
    [2024-03-27 15:58:12] TASK: Traversing tree to determine classification method.
    [2024-03-27 15:58:12] INFO: Completed 30 genomes in 0.02 seconds (1,364.18 genomes/second).
    [2024-03-27 15:58:15] TASK: Calculating average nucleotide identity using skani (v0.2.1).
    [2024-03-27 15:58:16] INFO: Completed 505 comparisons in 1.32 seconds (381.38 comparisons/second).
    [2024-03-27 15:58:19] INFO: 0 genome(s) have been classified using skani and pplacer.
    [2024-03-27 15:58:19] TASK: Placing 27 bacterial genomes into class-level reference tree 2 (2/8) with pplacer using 64 CPUs (be patient).
    [2024-03-27 16:05:53] INFO: Calculating RED values based on reference tree.
    [2024-03-27 16:05:56] TASK: Traversing tree to determine classification method.
    [2024-03-27 16:05:56] INFO: Completed 27 genomes in 0.04 seconds (606.99 genomes/second).
    [2024-03-27 16:05:59] TASK: Calculating average nucleotide identity using skani (v0.2.1).
    [2024-03-27 16:06:00] INFO: Completed 317 comparisons in 1.07 seconds (297.29 comparisons/second).
    [2024-03-27 16:06:03] INFO: 0 genome(s) have been classified using skani and pplacer.
    [2024-03-27 16:06:03] TASK: Placing 26 bacterial genomes into class-level reference tree 6 (3/8) with pplacer using 64 CPUs (be patient).
    [2024-03-27 16:12:28] INFO: Calculating RED values based on reference tree.
    [2024-03-27 16:12:30] TASK: Traversing tree to determine classification method.
    [2024-03-27 16:12:30] INFO: Completed 26 genomes in 0.05 seconds (497.16 genomes/second).
    [2024-03-27 16:12:31] TASK: Calculating average nucleotide identity using skani (v0.2.1).
    [2024-03-27 16:12:32] INFO: Completed 41 comparisons in 0.83 seconds (49.26 comparisons/second).
    [2024-03-27 16:12:33] INFO: 0 genome(s) have been classified using skani and pplacer.
    [2024-03-27 16:12:34] TASK: Placing 22 bacterial genomes into class-level reference tree 7 (4/8) with pplacer using 64 CPUs (be patient).
    [2024-03-27 16:18:24] INFO: Calculating RED values based on reference tree.
    [2024-03-27 16:18:27] TASK: Traversing tree to determine classification method.
    [2024-03-27 16:18:27] INFO: Completed 22 genomes in 0.03 seconds (715.55 genomes/second).
    [2024-03-27 16:18:28] TASK: Calculating average nucleotide identity using skani (v0.2.1).
    [2024-03-27 16:18:28] INFO: Completed 117 comparisons in 0.84 seconds (138.63 comparisons/second).
    [2024-03-27 16:18:30] INFO: 0 genome(s) have been classified using skani and pplacer.
    [2024-03-27 16:18:30] TASK: Placing 22 bacterial genomes into class-level reference tree 1 (5/8) with pplacer using 64 CPUs (be patient).
    [2024-03-27 16:26:01] INFO: Calculating RED values based on reference tree.
    [2024-03-27 16:26:04] TASK: Traversing tree to determine classification method.
    [2024-03-27 16:26:04] INFO: Completed 22 genomes in 0.05 seconds (486.20 genomes/second).
    [2024-03-27 16:26:05] TASK: Calculating average nucleotide identity using skani (v0.2.1).
    [2024-03-27 16:26:06] INFO: Completed 373 comparisons in 1.05 seconds (354.70 comparisons/second).
    [2024-03-27 16:26:08] INFO: 0 genome(s) have been classified using skani and pplacer.
    [2024-03-27 16:26:08] TASK: Placing 8 bacterial genomes into class-level reference tree 4 (6/8) with pplacer using 64 CPUs (be patient).
    [2024-03-27 16:32:45] INFO: Calculating RED values based on reference tree.
    [2024-03-27 16:32:48] TASK: Traversing tree to determine classification method.
    [2024-03-27 16:32:48] INFO: Completed 8 genomes in 0.15 seconds (52.49 genomes/second).
    [2024-03-27 16:32:48] TASK: Calculating average nucleotide identity using skani (v0.2.1).
    [2024-03-27 16:32:49] INFO: Completed 176 comparisons in 0.90 seconds (195.38 comparisons/second).
    [2024-03-27 16:32:50] INFO: 0 genome(s) have been classified using skani and pplacer.
    [2024-03-27 16:32:50] TASK: Placing 4 bacterial genomes into class-level reference tree 8 (7/8) with pplacer using 64 CPUs (be patient).
    [2024-03-27 16:35:30] INFO: Calculating RED values based on reference tree.
    [2024-03-27 16:35:31] TASK: Traversing tree to determine classification method.
    [2024-03-27 16:35:31] INFO: Completed 4 genomes in 0.00 seconds (5,959.93 genomes/second).
    [2024-03-27 16:35:32] TASK: Calculating average nucleotide identity using skani (v0.2.1).
    [2024-03-27 16:35:33] INFO: Completed 31 comparisons in 0.93 seconds (33.24 comparisons/second).
    [2024-03-27 16:35:33] INFO: 0 genome(s) have been classified using skani and pplacer.
    [2024-03-27 16:35:33] TASK: Placing 4 bacterial genomes into class-level reference tree 5 (8/8) with pplacer using 64 CPUs (be patient).
    [2024-03-27 16:40:57] INFO: Calculating RED values based on reference tree.
    [2024-03-27 16:40:59] TASK: Traversing tree to determine classification method.
    [2024-03-27 16:40:59] INFO: Completed 4 genomes in 0.00 seconds (4,607.86 genomes/second).
    [2024-03-27 16:40:59] TASK: Calculating average nucleotide identity using skani (v0.2.1).
    [2024-03-27 16:41:00] INFO: Completed 46 comparisons in 0.86 seconds (53.34 comparisons/second).
    [2024-03-27 16:41:00] INFO: 0 genome(s) have been classified using skani and pplacer.
    [2024-03-27 16:41:00] WARNING: 5 of 500 genomes have a warning (see summary file).
    [2024-03-27 16:41:00] INFO: Note that Tk classification mode is insufficient for publication of new taxonomic designations. New designations should be based on one or more de novo trees, an example of which can be produced by Tk in de novo mode.
    [2024-03-27 16:41:00] INFO: Done.