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
    * intermediate_results
        * :ref:`[prefix].[domain].classification_pplacer.tsv <files/classification_pplacer.tsv>`
        * :ref:`[prefix].[domain].classify.tree <files/classify.tree>`
        * pplacer
            * :ref:`pplacer.[domain].json <files/pplacer.domain.json>`
            * :ref:`pplacer.[domain].out <files/pplacer.domain.out>`
            * :ref:`[prefix].[domain].red_dictionary.tsv <files/red_dictionary.tsv>`
* :ref:`[prefix].[domain].summary.tsv <files/summary.tsv>`
* :ref:`[prefix].log <files/gtdbtk.log>`
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

    [2022-04-11 12:02:06] INFO: GTDB-Tk v2.0.0
    [2022-04-11 12:02:06] INFO: gtdbtk classify --genome_dir /tmp/gtdbtk/genomes --align_dir /tmp/gtdbtk/align --out_dir /tmp/gtdbtk/classify -x gz --cpus 2
    [2022-04-11 12:02:06] INFO: Using GTDB-Tk reference data version r207: /srv/db/gtdbtk/official/release207
    [2022-04-11 12:02:07] TASK: Placing 2 archaeal genomes into reference tree with pplacer using 2 CPUs (be patient).
    [2022-04-11 12:02:07] INFO: pplacer version: v1.1.alpha19-0-g807f6f3
    [2022-04-11 12:07:06] INFO: Calculating RED values based on reference tree.
    [2022-04-11 12:07:06] TASK: Traversing tree to determine classification method.
    [2022-04-11 12:07:06] INFO: Completed 2 genomes in 0.00 seconds (18,558.87 genomes/second).
    [2022-04-11 12:07:06] TASK: Calculating average nucleotide identity using FastANI (v1.32).
    [2022-04-11 12:07:08] INFO: Completed 4 comparisons in 1.61 seconds (2.49 comparisons/second).
    [2022-04-11 12:07:08] INFO: 2 genome(s) have been classified using FastANI and pplacer.
    [2022-04-11 12:07:08] INFO: Note that Tk classification mode is insufficient for publication of new taxonomic designations. New designations should be based on one or more de novo trees, an example of which can be produced by Tk in de novo mode.
    [2022-04-11 12:07:08] INFO: Done.
