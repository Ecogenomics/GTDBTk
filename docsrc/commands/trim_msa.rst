.. _commands/trim_msa:

trim_msa
========

The `trim_msa` command will trim a MSA given a user-specified mask file, or the archaeal/bacterial 
mask present in the reference data.

Arguments
---------

.. argparse::
   :module: gtdbtk.argparse
   :func: get_main_parser
   :prog: gtdbtk
   :path: trim_msa
   :nodefaultconst:


Example
-------
 

Input
^^^^^


.. code-block:: bash

    gtdbtk trim_msa --untrimmed_msa msa.faa --output msa_trim.faa --mask_file mask.txt



#### msa.faa

.. code-block:: text

    >genome_a
    AKLAK



#### mask.txt

.. code-block:: text
    
    01011


Output
^^^^^^


.. code-block:: text
    
    [2020-04-13 10:25:13] INFO: GTDB-Tk v1.1.0
    [2020-04-13 10:25:13] INFO: gtdbtk trim_msa --untrimmed_msa msa.faa --output msa_trim.faa --mask_file mask.txt
    [2020-04-13 10:25:13] INFO: Using GTDB-Tk reference data version r89: /release89
    [2020-04-13 10:25:13] INFO: Done.



#### msa_trim.faa

.. code-block:: text

    >genome_a
    KAK



