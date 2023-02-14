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

     gtdbtk classify --genome_dir genomes/ --align_dir align_output/ --out_dir classify_output --cpus 3



Output
^^^^^^


.. code-block:: text

    [2023-02-08 12:53:42] INFO: GTDB-Tk v2.2.0
    [2023-02-08 12:53:42] INFO: gtdbtk classify --align_dir align_3lines/ --batchfile 3lines_batchfile.tsv --out_dir 3classify_ani --mash_db mash_db_dir/ --cpus 20
    [2023-02-08 12:53:42] INFO: Using GTDB-Tk reference data version r207: /path/to/gtdbtk/database/release207_v2/
    [2023-02-08 12:53:43] INFO: Loading reference genomes.
    [2023-02-08 12:53:43] INFO: Using Mash version 2.2.2
    [2023-02-08 12:53:43] INFO: Loading data from existing Mash sketch file: 3classify_ani/classify/ani_screen/intermediate_results/mash/gtdbtk.user_query_sketch.msh
    [2023-02-08 12:53:43] INFO: Loading data from existing Mash sketch file: mash_db_dir/gtdb_ref_sketch.msh
    [2023-02-08 12:53:46] INFO: Calculating Mash distances.
    [2023-02-08 12:53:49] INFO: Calculating ANI with FastANI v1.3.
    [2023-02-08 12:53:49] INFO: Completed 12 comparisons in 0.44 seconds (27.54 comparisons/second).
    [2023-02-08 12:53:49] INFO: 2 genome(s) have been classified using the ANI pre-screening step.
    [2023-02-08 12:53:49] TASK: Placing 1 bacterial genomes into backbone reference tree with pplacer using 20 CPUs (be patient).
    [2023-02-08 12:53:49] INFO: pplacer version: v1.1.alpha19-0-g807f6f3
    [2023-02-08 12:55:02] INFO: Calculating RED values based on reference tree.
    [2023-02-08 12:55:03] INFO: 1 out of 1 have an class assignments. Those genomes will be reclassified.
    [2023-02-08 12:55:03] TASK: Placing 1 bacterial genomes into class-level reference tree 5 (1/1) with pplacer using 20 CPUs (be patient).
    [2023-02-08 12:57:38] INFO: Calculating RED values based on reference tree.
    [2023-02-08 12:57:40] TASK: Traversing tree to determine classification method.
    [2023-02-08 12:57:40] INFO: Completed 1 genome in 0.04 seconds (23.86 genomes/second).
    [2023-02-08 12:57:40] INFO: 0 genome(s) have been classified using FastANI and pplacer.
    [2023-02-08 12:57:40] WARNING: 1 of 3 genome has a warning (see summary file).
    [2023-02-08 12:57:40] INFO: Note that Tk classification mode is insufficient for publication of new taxonomic designations. New designations should be based on one or more de novo trees, an example of which can be produced by Tk in de novo mode.
    [2023-02-08 12:57:40] INFO: Done.