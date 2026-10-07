.. _commands/remove_labels:

remove_labels
=============

Remove the internal node labels (support values and taxon labels) from a Newick tree, to improve compatibility with
tree viewers.

Arguments
---------

.. argparse::
   :module: gtdbtk.cli
   :func: get_main_parser
   :prog: gtdbtk
   :path: remove_labels
   :nodefaultconst:


Example
-------

.. code-block:: bash

    gtdbtk remove_labels --input_tree gtdbtk.bac120.classify.tree --output_tree no_labels.tree
