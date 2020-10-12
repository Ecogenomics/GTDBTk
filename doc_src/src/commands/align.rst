.. _commands/align:

align
=====

Create a multiple sequence alignment based on the AR122/BAC120 marker set.


Arguments
---------

.. argparse::
   :module: gtdbtk.argparse
   :func: get_main_parser
   :prog: gtdbtk
   :path: align
   :nodefaultconst:


Files output
------------

* :ref:`[prefix].[domain].filtered.tsv <files/filtered.tsv>`
* :ref:`[prefix].[domain].msa.fasta <files/msa.fasta>`
* :ref:`[prefix].[domain].user_msa.fasta <files/user_msa.fasta>`
* :ref:`[prefix].log <files/gtdbtk.log>`
* :ref:`[prefix].warnings.log <files/gtdbtk.warnings.log>`
* :ref:`align/intermediate_results/[prefix].[domain].marker_info.tsv <files/marker_info.tsv>`


Example
-------


Input
^^^^^

.. code-block:: bash

    gtdbtk align --identify_dir identify_output/ --out_dir align_output --cpus 3


Output
^^^^^^

.. code-block:: text

    [2020-04-14 09:14:44] INFO: GTDB-Tk v1.1.0
    [2020-04-14 09:14:44] INFO: gtdbtk align --identify_dir identify_output/ --out_dir align_output --cpus 3
    [2020-04-14 09:14:44] INFO: Using GTDB-Tk reference data version r89: /release89
    [2020-04-14 09:14:44] INFO: Aligning markers in 3 genomes with 3 threads.
    [2020-04-14 09:14:44] INFO: Processing 3 genomes identified as archaeal.
    [2020-04-14 09:14:44] INFO: Read concatenated alignment for 1248 GTDB genomes.
    ==> Finished aligning 3 of 3 (100.0%) genomes.
    [2020-04-14 09:14:49] INFO: Masking columns of multiple sequence alignment using canonical mask.
    [2020-04-14 09:14:52] INFO: Masked alignment from 32675 to 5124 AAs.
    [2020-04-14 09:14:52] INFO: 0 user genomes have amino acids in <10.0% of columns in filtered MSA.
    [2020-04-14 09:14:52] INFO: Creating concatenated alignment for 1251 GTDB and user genomes.
    [2020-04-14 09:14:52] INFO: Creating concatenated alignment for 3 user genomes.
    [2020-04-14 09:14:52] INFO: Done.
