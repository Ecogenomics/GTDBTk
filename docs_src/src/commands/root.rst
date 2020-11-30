.. _commands/root:

root
====

Root a tree using an outgroup.

Arguments
---------

.. argparse::
   :module: gtdbtk.cli
   :func: get_main_parser
   :prog: gtdbtk
   :path: root
   :nodefaultconst:


Example
-------


Input

    
.. code-block:: bash

    gtdbtk root --input_tree input.tree --outgroup_taxon p__Nanoarchaeota --output_tree output.tree
    


Output
^^^^^^


.. code-block:: text

    [2020-04-14 08:26:53] INFO: GTDB-Tk v1.1.0
    [2020-04-14 08:26:53] INFO: gtdbtk root --input_tree input.tree --outgroup_taxon p__Nanoarchaeota --output_tree output.tree
    [2020-04-14 08:26:53] INFO: Using GTDB-Tk reference data version r89: /release89
    [2020-04-14 08:26:53] INFO: Identifying genomes from the specified outgroup.
    [2020-04-14 08:26:53] INFO: Identified 101 outgroup taxa in the tree.
    [2020-04-14 08:26:53] INFO: Identified 1151 ingroup taxa in the tree.
    [2020-04-14 08:26:53] INFO: Outgroup is monophyletic.
    [2020-04-14 08:26:53] INFO: Rerooting tree.
    [2020-04-14 08:26:53] INFO: Rerooted tree written to: output.tree
    [2020-04-14 08:26:53] INFO: Done.
    
