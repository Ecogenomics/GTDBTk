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



Example
-------

Input
^^^^^

.. code-block:: bash

    gtdbtk ani_rep --batchfile genomes/500_batchfile.tsv -x fa --out_dir test_ani_reps --cpus 90


Output
^^^^^^

.. code-block:: text

    [2025-08-05 17:13:33] INFO: GTDB-Tk v2.5.0
    [2025-08-05 17:13:33] INFO: gtdbtk ani_rep --batchfile genomes/500_batchfile.tsv -x fa --out_dir test_ani_reps --cpus 90
    [2025-08-05 17:13:33] INFO: Using GTDB-Tk reference data version r226: /srv/db/gtdbtk/official/release226
    [2025-08-05 17:13:34] INFO: Loading reference genomes.
    [2025-08-05 17:13:39] INFO: Calculating all vs all ANI with skani v0.2.1.
    [2025-08-05 17:13:40] INFO: Sketching genomes
    [2025-08-05 17:15:52] INFO: Sketches done: 2min 11secs
    [2025-08-05 17:15:52] INFO: Running comparisons
    [2025-08-05 17:21:15] INFO: Comparisons finished, capturing results.
    [2025-08-05 17:23:26] INFO: Summary of results saved to: test_ani_reps/gtdbtk.ani_summary.tsv
    [2025-08-05 17:23:27] INFO: Closest representative hits saved to: test_ani_reps/gtdbtk.ani_closest.tsv
    [2025-08-05 17:23:27] INFO: Done.



