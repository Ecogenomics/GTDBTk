.. _commands/classify:

classify
========

Determine taxonomic classification of genomes.

Arguments
---------

.. argparse::
   :module: gtdbtk.argparse
   :func: get_main_parser
   :prog: gtdbtk
   :path: classify
   :nodefaultconst:


Files output
------------

* classify
    * intermediate_results
        * :ref:`[prefix].[domain].classification_pplacer.tsv <files/classification_pplacer.tsv>`
        * :ref:`[prefix].[domain].red_dictionary.tsv <files/red_dictionary.tsv>`
        * pplacer
            * :ref:`pplacer.[domain].json <files/pplacer.domain.json>`
            * :ref:`pplacer.[domain].out <files/pplacer.domain.out>`
* :ref:`[prefix].[domain].classify.tree <files/classify.tree>`
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

    [2020-04-14 09:22:35] INFO: GTDB-Tk v1.1.0
    [2020-04-14 09:22:35] INFO: gtdbtk classify --genome_dir genomes/ --align_dir align_output/ --out_dir classify_output --cpus 3
    [2020-04-14 09:22:35] INFO: Using GTDB-Tk reference data version r89: /release89
    [2020-04-14 09:22:35] INFO: Placing 3 archaeal genomes into reference tree with pplacer using 3 cpus (be patient).
    Placing genomes |##################################################| 3/3 (100.00%)
    [2020-04-14 09:23:32] INFO: pplacer version: v1.1.alpha19-0-g807f6f3
    [2020-04-14 09:23:33] INFO: Calculating average nucleotide identity using FastANI.
    [2020-04-14 09:23:33] INFO: fastANI version: 1.3
    ==> Processing 24 of 24 (100.0%) comparisons.
    [2020-04-14 09:23:38] INFO: 3 genome(s) have been classified using FastANI and pplacer.
    [2020-04-14 09:23:38] INFO: Calculating RED values based on reference tree.
    [2020-04-14 09:23:38] INFO: Done.
