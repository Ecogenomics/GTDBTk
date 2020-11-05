.. _commands/check_install:

check_install
=============

The `check_install` command is used to verify the integrity of the GTDB-Tk reference data.

If any inconsistencies are identified then the program will exit with code 1 and 
`HASH MISMATCH` will be displayed next to the inconsistent item.

Arguments
---------

.. argparse::
   :module: gtdbtk.argparse
   :func: get_main_parser
   :prog: gtdbtk
   :path: check_install
   :nodefaultconst:


Example
-------

Input
^^^^^

.. code-block:: bash

    gtdbtk check_install


Output
^^^^^^

.. code-block:: text
    
    [2020-11-04 09:35:16] INFO: GTDB-Tk v1.4.0
    [2020-11-04 09:35:16] INFO: gtdbtk check_install
    [2020-11-04 09:35:16] INFO: Using GTDB-Tk reference data version r95: /release95
    [2020-11-04 09:35:16] INFO: Running install verification
    [2020-11-04 09:35:16] INFO: Checking that all third-party software are on the system path:
    [2020-11-04 09:35:16] INFO:          |-- FastTree         OK
    [2020-11-04 09:35:16] INFO:          |-- FastTreeMP       OK
    [2020-11-04 09:35:16] INFO:          |-- fastANI          OK
    [2020-11-04 09:35:16] INFO:          |-- guppy            OK
    [2020-11-04 09:35:16] INFO:          |-- hmmalign         OK
    [2020-11-04 09:35:16] INFO:          |-- hmmsearch        OK
    [2020-11-04 09:35:16] INFO:          |-- mash             OK
    [2020-11-04 09:35:16] INFO:          |-- pplacer          OK
    [2020-11-04 09:35:16] INFO:          |-- prodigal         OK
    [2020-11-04 09:35:16] INFO: Checking /release95
    [2020-11-04 09:35:16] INFO:          |-- pplacer          OK
    [2020-11-04 09:35:16] INFO:          |-- masks            OK
    [2020-11-04 09:35:17] INFO:          |-- markers          OK
    [2020-11-04 09:35:17] INFO:          |-- radii            OK
    [2020-11-04 09:35:20] INFO:          |-- msa              OK
    [2020-11-04 09:35:20] INFO:          |-- metadata         OK
    [2020-11-04 09:35:20] INFO:          |-- taxonomy         OK
    [2020-11-04 09:47:36] INFO:          |-- fastani          OK
    [2020-11-04 09:47:36] INFO:          |-- mrca_red         OK
    [2020-11-04 09:47:36] INFO: Done.
