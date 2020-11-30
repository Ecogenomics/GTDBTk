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

    [2020-04-13 10:51:58] INFO: GTDB-Tk v1.1.0
    [2020-04-13 10:51:58] INFO: gtdbtk ani_rep --genome_dir genomes/ --out_dir ani_rep/ --cpus 70
    [2020-04-13 10:51:58] INFO: Using GTDB-Tk reference data version r89: /release89
    [2020-04-13 10:51:59] INFO: Using Mash version 2.2.2
    [2020-04-13 10:51:59] INFO: Creating Mash sketch file: ani_rep/intermediate_results/mash/gtdbtk.user_query_sketch.msh
    ==> Sketching 3 of 3 (100.0%) genomes
    [2020-04-13 10:51:59] INFO: Creating Mash sketch file: ani_rep/intermediate_results/mash/gtdbtk.gtdb_ref_sketch.msh
    ==> Sketching 24706 of 24706 (100.0%) genomes
    [2020-04-13 10:53:13] INFO: Calculating Mash distances.
    [2020-04-13 10:53:14] INFO: Calculating ANI with FastANI.
    ==> Processing 874 of 874 (100.0%) comparisons.
    [2020-04-13 10:53:23] INFO: Summary of results saved to: ani_rep/gtdbtk.ani_summary.tsv
    [2020-04-13 10:53:23] INFO: Closest representative hits saved to: ani_rep/gtdbtk.ani_closest.tsv
    [2020-04-13 10:53:23] INFO: Done.



