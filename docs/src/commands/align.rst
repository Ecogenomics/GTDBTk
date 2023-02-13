.. _commands/align:

align
=====

Create a multiple sequence alignment based on the AR53/BAC120 marker set.


Arguments
---------

.. argparse::
   :module: gtdbtk.cli
   :func: get_main_parser
   :prog: gtdbtk
   :path: align
   :nodefaultconst:


Files output
------------


* :ref:`[prefix].log <files/gtdbtk.log>`
* :ref:`[prefix].json <files/gtdbtk.json>`
* :ref:`[prefix].warnings.log <files/gtdbtk.warnings.log>`
* align
    * :ref:`[prefix].[domain].msa.fasta.gz <files/msa.fasta>`
    * :ref:`[prefix].[domain].user_msa.fasta.gz <files/user_msa.fasta>`
    * :ref:`[prefix].[domain].filtered.tsv <files/filtered.tsv>`
    * intermediate_results
        * :ref:`[prefix].[domain].marker_info.tsv <files/marker_info.tsv>`

Example
-------


Input
^^^^^

.. code-block:: bash

    gtdbtk align --identify_dir identify_output/ --out_dir align_output --cpus 3


Output
^^^^^^

.. code-block:: text

    [2022-04-11 11:59:14] INFO: GTDB-Tk v2.0.0
    [2022-04-11 11:59:14] INFO: gtdbtk align --identify_dir /tmp/gtdbtk/identify --out_dir /tmp/gtdbtk/align --cpus 2
    [2022-04-11 11:59:14] INFO: Using GTDB-Tk reference data version r207: /srv/db/gtdbtk/official/release207
    [2022-04-11 11:59:15] INFO: Aligning markers in 3 genomes with 2 CPUs.
    [2022-04-11 11:59:16] INFO: Processing 3 genomes identified as archaeal.
    [2022-04-11 11:59:16] INFO: Read concatenated alignment for 3,412 GTDB genomes.
    [2022-04-11 11:59:16] TASK: Generating concatenated alignment for each marker.
    [2022-04-11 11:59:16] INFO: Completed 3 genomes in 0.01 seconds (139.73 genomes/second).
    [2022-04-11 11:59:16] TASK: Aligning 52 identified markers using hmmalign 3.1b2 (February 2015).
    [2022-04-11 11:59:17] INFO: Completed 52 markers in 0.86 seconds (60.66 markers/second).
    [2022-04-11 11:59:17] TASK: Masking columns of archaeal multiple sequence alignment using canonical mask.
    [2022-04-11 11:59:21] INFO: Completed 3,414 sequences in 4.19 seconds (815.22 sequences/second).
    [2022-04-11 11:59:21] INFO: Masked archaeal alignment from 13,540 to 10,153 AAs.
    [2022-04-11 11:59:21] INFO: 0 archaeal user genomes have amino acids in <10.0% of columns in filtered MSA.
    [2022-04-11 11:59:21] INFO: Creating concatenated alignment for 3,414 archaeal GTDB and user genomes.
    [2022-04-11 11:59:23] INFO: Creating concatenated alignment for 3 archaeal user genomes.
    [2022-04-11 11:59:23] INFO: Done.
