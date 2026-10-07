.. _files/failed_genomes.tsv:

failed_genomes.tsv
==================

Genomes excluded from the analysis because Prodigal failed to call any genes, or because the genome file was empty.
Tab-separated: genome identifier and reason. These genomes are reported as ``Unclassified`` in the bac120
:ref:`summary file <files/summary.tsv>`.

Produced by
-----------

* :ref:`commands/identify`
* :ref:`commands/classify_wf`
* :ref:`commands/de_novo_wf`

Example
-------

.. code-block:: text

    GCA_000002165.1	No genes were called by Prodigal
    GCA_000002175.1	No genes were called by Prodigal
    GCA_000002185.1	Empty file
