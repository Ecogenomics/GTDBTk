.. _commands/convert_to_itol:

convert_to_itol
===============

The `convert_to_itol` command will remove internal labels from Newick tree, making it suitable for visualization in `iTOL <http://itol.embl.de/>`_.  

Arguments
---------

.. argparse::
   :module: gtdbtk.cli
   :func: get_main_parser
   :prog: gtdbtk
   :path: convert_to_itol
   :nodefaultconst:

Example
-------

Input
^^^^^


.. code-block:: bash

    gtdbtk convert_to_itol --input some_tree.tree --output itol.tree


Output
^^^^^^


.. code-block:: text

    [2022-06-30 18:44:54] INFO: GTDB-Tk v2.1.0
    [2022-06-30 18:44:54] INFO: gtdbtk convert_to_itol --input /tmp/decorated.tree --output new.tree
    [2022-06-30 18:44:54] INFO: Using GTDB-Tk reference data version r207: /gtdbtk-data
    [2022-06-30 18:44:54] INFO: Convert GTDB-Tk tree to iTOL format
    [2022-06-30 18:44:54] INFO: Done.

