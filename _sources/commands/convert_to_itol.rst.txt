.. _commands/convert_to_itol:

convert_to_itol
===============

The ``convert_to_itol`` command reformats a GTDB-Tk Newick tree for visualization in `iTOL <https://itol.embl.de/>`_:
taxon labels on internal nodes are kept (with ``;`` replaced by ``|``), and support values are moved into square
brackets after the branch lengths. To remove all internal labels instead, use :ref:`commands/remove_labels`.

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

    gtdbtk convert_to_itol --input_tree some_tree.tree --output_tree itol.tree


Output
^^^^^^


.. code-block:: text

    [2022-06-30 18:44:54] INFO: GTDB-Tk v2.1.0
    [2022-06-30 18:44:54] INFO: gtdbtk convert_to_itol --input_tree /tmp/decorated.tree --output_tree new.tree
    [2022-06-30 18:44:54] INFO: Using GTDB-Tk reference data version r207: /gtdbtk-data
    [2022-06-30 18:44:54] INFO: Convert GTDB-Tk tree to iTOL format
    [2022-06-30 18:44:54] INFO: Done.

