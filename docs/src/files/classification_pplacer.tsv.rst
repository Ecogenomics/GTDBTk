.. _files/classification_pplacer.tsv:

classification_pplacer.tsv
==========================

The classification of query genomes based only on their placement in the reference tree.

* ``[prefix].ar53.classification_pplacer.tsv``, ``[prefix].bac120.class_level.classification_pplacer_tree_[index].tsv``
  and, with ``--full_tree``, ``[prefix].bac120.classification_pplacer.tsv``: two tab-separated columns without a header,
  the genome identifier and the taxonomy from its placement (first example below).
* ``[prefix].bac120.backbone.classification_pplacer.tsv`` (split mode): placement on the bacterial backbone tree, with a
  header and the columns ``user_genome``, ``gtdb_taxonomy_red``, ``gtdb_taxonomy_terminal``, ``pplacer_taxonomy``,
  ``is_terminal`` and ``red``.

Produced by
-----------

 * :ref:`commands/classify`
 * :ref:`commands/classify_wf`


Example
-------

.. code-block:: text

    genome_2	d__Archaea;p__Thermoplasmatota;c__Thermoplasmata;o__Methanomassiliicoccales;f__Methanomethylophilaceae;g__VadinCA11;s__
    genome_3	d__Archaea;p__Thermoplasmatota;c__Thermoplasmata;o__Methanomassiliicoccales;f__Methanomethylophilaceae;g__VadinCA11;s__
    genome_1	d__Archaea;p__Methanobacteriota;c__Methanobacteria;o__Methanobacteriales;f__Methanobacteriaceae;g__Methanobrevibacter;s__
