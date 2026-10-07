.. _files/ani_summary.tsv:

ani_summary.tsv
===============

ANI and alignment fraction (AF) between the user genomes and the GTDB species representative genomes, computed with
skani. Tab-separated with the following columns:

* ``user_genome``: identifier of the query genome.
* ``reference_genome``: accession of the GTDB species representative genome.
* ``skani_ani``: ANI between the query and the reference genome.
* ``skani_af``: AF between the query and the reference genome (the larger of the query and reference alignment
  fractions reported by skani).
* ``reference_taxonomy``: GTDB taxonomy of the reference genome.
* ``other_related_references(genome_id,species_name,radius,ANI,AF)``: ``classify_wf`` only, see
  :ref:`summary.tsv <files/summary.tsv>`.

Two versions of this file are written:

* :ref:`commands/ani_rep` writes ``<prefix>.ani_summary.tsv`` with one row for every skani hit of every query genome
  (sorted by decreasing AF, then ANI).
* :ref:`commands/classify_wf` writes ``classify/ani_screen/<prefix>.<domain>.ani_summary.tsv`` with one row per genome
  assigned to a species by the ANI screen (the representative it was assigned to), plus the
  ``other_related_references`` column. It is only written if at least one genome was assigned.

Produced by
-----------

* :ref:`commands/ani_rep`
* :ref:`commands/classify_wf`


Example
-------

.. code-block:: text

    user_genome	reference_genome	skani_ani	skani_af	reference_taxonomy
    genome_1	GCF_000024185.1	100.0	1.0	d__Archaea;p__Methanobacteriota;c__Methanobacteria;o__Methanobacteriales;f__Methanobacteriaceae;g__Methanobrevibacter;s__Methanobrevibacter ruminantium
    genome_1	GCA_900321995.1	80.9	0.7	d__Archaea;p__Methanobacteriota;c__Methanobacteria;o__Methanobacteriales;f__Methanobacteriaceae;g__Methanobrevibacter;s__Methanobrevibacter sp900321995
    genome_1	GCF_900114585.1	79.96	0.55	d__Archaea;p__Methanobacteriota;c__Methanobacteria;o__Methanobacteriales;f__Methanobacteriaceae;g__Methanobrevibacter;s__Methanobrevibacter olleyae
    genome_2	GCA_002498365.1	99.16	0.94	d__Archaea;p__Thermoplasmatota;c__Thermoplasmata;o__Methanomassiliicoccales;f__Methanomethylophilaceae;g__VadinCA11;s__VadinCA11 sp002498365
    genome_2	GCA_002505345.1	89.92	0.89	d__Archaea;p__Thermoplasmatota;c__Thermoplasmata;o__Methanomassiliicoccales;f__Methanomethylophilaceae;g__VadinCA11;s__VadinCA11 sp002505345
