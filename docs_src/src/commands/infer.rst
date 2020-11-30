.. _commands/infer:

infer
=====

Infer tree from multiple sequence alignment.

Arguments
---------

.. argparse::
   :module: gtdbtk.cli
   :func: get_main_parser
   :prog: gtdbtk
   :path: infer
   :nodefaultconst:


## Files output

* :ref:`[prefix].log <files/gtdbtk.log>`
* :ref:`[prefix].unrooted.tree <files/unrooted.tree>`
* :ref:`[prefix].warnings.log <files/gtdbtk.warnings.log>`
* infer/intermediate_results/
    * :ref:`[prefix].fasttree.log <files/fasttree.log>`
    * :ref:`[prefix].tree.log <files/tree.log>`


Example
-------


Input
^^^^^


.. code-block:: bash

    gtdbtk infer --msa_file msa.faa --out_dir infer_out



Output
^^^^^^

.. code-block:: text
    
    [2020-04-14 09:37:55] INFO: GTDB-Tk v1.1.0
    [2020-04-14 09:37:55] INFO: gtdbtk infer --msa_file msa.faa --out_dir infer_out
    [2020-04-14 09:37:55] INFO: Using GTDB-Tk reference data version r89: /release89
    [2020-04-14 09:37:55] INFO: Inferring FastTree (WAG, +gamma, support) using a maximum of 1 CPUs.
    [2020-04-14 09:37:55] INFO: FastTree version: 2.1.10
    [2020-04-14 09:37:55] INFO: FastTree version: 2.1.10
    [2020-04-14 09:37:55] INFO: Done.

