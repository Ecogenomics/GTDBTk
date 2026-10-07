.. _files/ani_closest.tsv:

ani_closest.tsv
===============

For each user genome, the closest GTDB species representative genome: the hit with the highest ANI among those with
an alignment fraction (AF) ≥ ``--min_af``. Tab-separated with the following columns:

* ``user_genome``: identifier of the query genome.
* ``reference_genome``: accession of the closest representative genome.
* ``skani_ani``: ANI between the query and the reference genome.
* ``skani_af``: AF between the query and the reference genome.
* ``reference_taxonomy``: GTDB taxonomy of the reference genome.
* ``satisfies_gtdb_circumscription_criteria``: ``True`` if the ANI is within the species-specific ANI circumscription
  radius of the reference genome and the AF is ≥ 0.5.

Genomes without any hit with AF ≥ ``--min_af`` are reported with ``no result``.

Produced by
-----------

* :ref:`commands/ani_rep`

Example
-------

.. code-block:: text

    user_genome	reference_genome	skani_ani	skani_af	reference_taxonomy	satisfies_gtdb_circumscription_criteria
    genome_1	GCF_000024185.1	100.0	1.0	d__Archaea;p__Methanobacteriota;c__Methanobacteria;o__Methanobacteriales;f__Methanobacteriaceae;g__Methanobrevibacter;s__Methanobrevibacter ruminantium	True
    genome_2	GCA_002498365.1	99.16	0.94	d__Archaea;p__Thermoplasmatota;c__Thermoplasmata;o__Methanomassiliicoccales;f__Methanomethylophilaceae;g__VadinCA11;s__VadinCA11 sp002498365	True
    genome_3	no result	no result	no result	no result	no result
