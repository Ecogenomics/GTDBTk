.. _commands/convert_to_species:

convert_to_species
==================

Replace the GTDB genome accessions at the leaves of a GTDB-Tk tree with their GTDB species names.

Leaves that are not GTDB genomes (e.g. user genomes) are left unchanged, unless they are listed in a
``--custom_taxonomy_file`` (same format as in :ref:`commands/de_novo_wf`); entries in that file take precedence over
the GTDB taxonomy. With ``--all_ranks``, the full 7-rank taxonomy is used as the leaf label instead of the species name.

Arguments
---------

.. argparse::
   :module: gtdbtk.cli
   :func: get_main_parser
   :prog: gtdbtk
   :path: convert_to_species
   :nodefaultconst:


Example
-------

.. code-block:: bash

    gtdbtk convert_to_species --input_tree gtdbtk.bac120.classify.tree --output_tree species.tree
