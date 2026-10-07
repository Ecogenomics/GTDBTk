.. _commands/infer_ranks:

infer_ranks
===========

Establish the taxonomic ranks of the internal nodes of a rooted tree using relative evolutionary divergence (RED).

The tree is scaled with RED from the ``--ingroup_taxon`` node, using the RED of that taxon in the GTDB reference tree.
Each internal node within the ingroup is then labelled with its RED value and the ranks whose median RED is within
0.1 of it, closest first, e.g. ``RED=0.762|family&genus``. The input tree must be rooted and the ingroup taxon must be labelled in it
(e.g. a tree produced by :ref:`commands/root` and :ref:`commands/decorate`).

Arguments
---------

.. argparse::
   :module: gtdbtk.cli
   :func: get_main_parser
   :prog: gtdbtk
   :path: infer_ranks
   :nodefaultconst:


Example
-------

.. code-block:: bash

    gtdbtk infer_ranks --input_tree decorated.tree --ingroup_taxon c__Bacilli --output_tree ranks.tree
