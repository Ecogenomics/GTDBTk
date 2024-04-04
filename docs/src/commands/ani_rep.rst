.. _commands/ani_rep:

ani_rep
=======

Compute the ANI of input genomes to all GTDB-Tk representative genomes.

Arguments
---------

.. argparse::
   :module: gtdbtk.cli
   :func: get_main_parser
   :prog: gtdbtk
   :path: ani_rep
   :nodefaultconst:


Files output
------------


* :ref:`[prefix].ani_closest.tsv <files/ani_closest.tsv>`
* :ref:`[prefix].ani_summary.tsv <files/ani_summary.tsv>`
* :ref:`[prefix].log <files/gtdbtk.log>`
* :ref:`[prefix].warnings.log <files/gtdbtk.warnings.log>`
* :ref:`[prefix].warnings.log <files/gtdbtk.warnings.log>`

* intermediate_results/mash/

    * :ref:`[prefix].gtdb_ref_sketch.msh <files/gtdbtk_ref_sketch.msh>`
    * :ref:`[prefix].mash_distances.msh <files/mash_distances.msh>`
    * :ref:`[prefix].user_query_sketch.msh <files/user_query_sketch.msh>`


Example
-------

Input
^^^^^

.. code-block:: bash

    gtdbtk ani_rep --genome_dir genomes/ --out_dir ani_rep/ --cpus 70


Output
^^^^^^

.. code-block:: text

    [2024-03-27 16:43:25] INFO: GTDB-Tk v2.3.2
    [2024-03-27 16:43:25] INFO: gtdbtk ani_rep --batchfile genomes/500_batchfile.tsv --out_dir user_vs_reps --cpus 90
    [2024-03-27 16:43:25] INFO: Using GTDB-Tk reference data version r214: /srv/db/gtdbtk/official/release214_skani/release214
    [2024-03-27 16:43:25] INFO: Loading reference genomes.
    [2024-03-27 16:43:25] INFO: Using Mash version 2.2.2
    [2024-03-27 16:43:25] INFO: Creating Mash sketch file: user_vs_reps/intermediate_results/mash/gtdbtk.user_query_sketch.msh
    [2024-03-27 16:43:27] INFO: Completed 500 genomes in 1.42 seconds (351.61 genomes/second).
    [2024-03-27 16:43:27] INFO: Creating Mash sketch file: user_vs_reps/intermediate_results/mash/gtdbtk.gtdb_ref_sketch.msh
    [2024-03-27 16:46:55] INFO: Completed 85,205 genomes in 3.47 minutes (24,519.48 genomes/minute).
    [2024-03-27 16:46:55] INFO: Calculating Mash distances.
    [2024-03-27 16:47:37] INFO: Calculating ANI with skani v0.2.1.
    [2024-03-27 16:47:45] INFO: Completed 4,383 comparisons in 7.68 seconds (570.58 comparisons/second).
    [2024-03-27 16:47:46] INFO: Summary of results saved to: user_vs_reps/gtdbtk.ani_summary.tsv
    [2024-03-27 16:47:46] INFO: Closest representative hits saved to: user_vs_reps/gtdbtk.ani_closest.tsv
    [2024-03-27 16:47:46] INFO: Done.



