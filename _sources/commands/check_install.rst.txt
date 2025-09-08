.. _commands/check_install:

check_install
=============

The `check_install` command is used to verify the integrity of the GTDB-Tk reference data.

If any inconsistencies are identified then the program will exit with code 1 and 
`HASH MISMATCH` will be displayed next to the inconsistent item.

Arguments
---------

.. argparse::
   :module: gtdbtk.cli
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
    
    [2025-08-05 17:02:59] INFO: GTDB-Tk v2.5.0
    [2025-08-05 17:02:59] INFO: gtdbtk check_install
    [2025-08-05 17:02:59] INFO: Using GTDB-Tk reference data version r226: /release226
    [2025-08-05 17:02:59] INFO: Running install verification
    [2025-08-05 17:02:59] INFO: Checking that all third-party software are on the system path:
    [2025-08-05 17:02:59] INFO:          |-- FastTree         OK
    [2025-08-05 17:02:59] INFO:          |-- FastTreeMP       OK
    [2025-08-05 17:02:59] INFO:          |-- guppy            OK
    [2025-08-05 17:02:59] INFO:          |-- hmmalign         OK
    [2025-08-05 17:02:59] INFO:          |-- hmmsearch        OK
    [2025-08-05 17:02:59] INFO:          |-- pplacer          OK
    [2025-08-05 17:02:59] INFO:          |-- prodigal         OK
    [2025-08-05 17:02:59] INFO:          |-- skani            OK
    [2025-08-05 17:02:59] INFO: Checking integrity of reference package: /release226
    [2025-08-05 17:02:59] INFO:          |-- pplacer          OK
    [2025-08-05 17:02:59] INFO:          |-- masks            OK
    [2025-08-05 17:03:00] INFO:          |-- markers          OK
    [2025-08-05 17:03:00] INFO:          |-- radii            OK
    [2025-08-05 17:03:05] INFO:          |-- msa              OK
    [2025-08-05 17:03:05] INFO:          |-- metadata         OK
    [2025-08-05 17:03:05] INFO:          |-- taxonomy         OK
    [2025-08-05 17:06:35] INFO:          |-- skani            OK
    [2025-08-05 17:06:35] INFO:          |-- mrca_red         OK
    [2025-08-05 17:06:35] INFO: Done.
