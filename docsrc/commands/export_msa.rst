.. _commands/export_msa:

export_msa
==========

The `export_msa` will export the untrimmed archaeal or bacterial MSA used in the reference data.

Arguments
---------

.. argparse::
   :module: gtdbtk.argparse
   :func: get_main_parser
   :prog: gtdbtk
   :path: export_msa
   :nodefaultconst:


Example
-------


Input
^^^^^

.. code-block:: bash

    gtdbtk export_msa --domain arc --output /tmp/msa.faa



Output
^^^^^^


.. code-block:: text
    
    [2020-04-13 10:03:05] INFO: GTDB-Tk v1.1.0
    [2020-04-13 10:03:05] INFO: gtdbtk export_msa --domain arc --output /tmp/msa.faa
    [2020-04-13 10:03:05] INFO: Using GTDB-Tk reference data version r89: /release89
    [2020-04-13 10:03:05] INFO: Done.
    
