.. _files/summary.tsv:


summary.tsv
===========

Classifications provided by GTDB-Tk are in the files ``<prefix>.bac120.summary.tsv`` and ``<prefix>.ar53.summary.tsv``
for bacterial and archaeal genomes, respectively. These are tab-separated files with one row per genome and the
following 20 columns. Empty values are written as ``N/A``.

Genomes are first compared with skani to all GTDB species representative genomes (ANI screen). A genome is assigned
to a species when the closest representative with an alignment fraction (AF) ≥ ``--min_af`` (default 0.5) has an ANI
within that representative's species-specific ANI circumscription radius. All other genomes are placed in the reference
tree with pplacer and classified from their placement and relative evolutionary divergence (RED). Genomes assigned by
the ANI screen are only placed in the tree when ``--place_species`` is used.

* ``user_genome``: unique identifier of the query genome, taken from the FASTA file name (or the batch file).
* ``classification``: GTDB taxonomy string inferred by GTDB-Tk. An unassigned species (``s__``) means that no
  reference genome with an AF ≥ ``--min_af`` has an ANI within its circumscription radius, or that the genome is placed
  outside a named genus. Genomes for which gene calling failed or no marker genes were found are reported as
  ``Unclassified`` (in the bac120 file), and genomes removed by ``--min_perc_aa`` as ``Unclassified Bacteria`` or
  ``Unclassified Archaea``, with the reason in ``warnings``.
* ``closest_genome_reference``: accession of the species representative genome to which the query genome was assigned
  based on ANI and AF. ``N/A`` when no species was assigned by ANI.
* ``closest_genome_reference_radius``: species-specific ANI circumscription radius of the above reference genome.
* ``closest_genome_taxonomy``: GTDB taxonomy of the above reference genome.
* ``closest_genome_ani``: ANI between the query and the above reference genome.
* ``closest_genome_af``: AF between the query and the above reference genome (the larger of the query and reference
  alignment fractions reported by skani).
* ``closest_placement_reference``: for genomes assigned by the ANI screen and placed in the tree with
  ``--place_species``, the accession of the reference genome on the terminal branch where the query genome is placed.
  ``N/A`` otherwise.
* ``closest_placement_radius``: species-specific ANI circumscription radius of the above reference genome.
* ``closest_placement_taxonomy``: GTDB taxonomy of the above reference genome.
* ``closest_placement_ani``: ANI between the query and the above reference genome, when available.
* ``closest_placement_af``: AF between the query and the above reference genome, when available.
* ``pplacer_taxonomy``: taxonomy inferred from the pplacer placement of the query genome in the reference tree.
* ``classification_method``: rule used to classify the genome. One of:

  * ``ani_screen``: species assigned from ANI and AF alone during the ANI screen; the genome is not placed in the
    reference tree (unless ``--place_species`` is used).
  * ``taxonomic classification defined by topology and ANI``: with ``--place_species``, a genome assigned by the ANI
    screen that was also placed in the reference tree; the species assignment comes from ANI and AF.
  * ``taxonomic classification fully defined by topology``: the classification follows directly from the placement of
    the genome in the reference tree.
  * ``taxonomic novelty determined using RED``: the placement and the RED value of the genome were used to determine
    the classification (e.g. a putative novel genus or family).

* ``note``: additional information about the classification. Several notes are separated by ``;``. Values:

  * ``classification based on ANI only``: species assigned by the ANI screen.
  * ``topological placement and ANI have congruent species assignments``: the placement reference and the closest
    reference by ANI/AF are the same species.
  * ``topological placement and ANI have incongruent species assignments``: they differ; the species assignment
    follows ANI/AF.
  * ``classification based on placement in backbone tree``: in split mode (bacteria, default), the genome was
    classified on the backbone tree only.
  * ``classification based on placement in class-level tree``: in split mode, the classification comes from the
    class-level tree.
  * ``classification based on consensus between backbone and class-level tree``: in split mode, the backbone and
    class-level placements were combined.

* ``other_related_references(genome_id,species_name,radius,ANI,AF)``: other reference genomes found close to the
  query genome by the ANI screen, separated by ``;``. Each entry gives the accession, species name, circumscription
  radius, ANI and AF. Only reported for genomes assigned to a species by ANI.
* ``msa_percent``: percentage of the multiple sequence alignment spanned by the genome (i.e. percentage of columns
  with an amino acid).
* ``translation_table``: translation table used by Prodigal to call genes (``11`` or ``4``, or the table given in the
  batch file).
* ``red_value``: relative evolutionary divergence (RED) of the query genome, when it was needed for the
  classification. Not calculated when the genome is classified by ANI.
* ``warnings``: unusual characteristics of the query genome that may affect the taxonomic assignment, separated by
  ``;`` (e.g. a high percentage of markers with multiple hits, a questionable domain, or a genome not assigned to the
  closest species because it falls outside that species' ANI circumscription radius).


Produced by
-----------

 * :ref:`commands/classify`
 * :ref:`commands/classify_wf`

Example
-------

``genome_2`` shows a genome assigned by the ANI screen and placed in the tree with ``--place_species``.

.. code-block:: text

    user_genome	classification	closest_genome_reference	closest_genome_reference_radius	closest_genome_taxonomy	closest_genome_ani	closest_genome_af	closest_placement_reference	closest_placement_radius	closest_placement_taxonomy	closest_placement_ani	closest_placement_af	pplacer_taxonomy	classification_method	note	other_related_references(genome_id,species_name,radius,ANI,AF)	msa_percent	translation_table	red_value	warnings
    genome_1	d__Archaea;p__Methanobacteriota;c__Methanobacteria;o__Methanobacteriales;f__Methanobacteriaceae;g__Methanobrevibacter;s__Methanobrevibacter ruminantium	GCF_000024185.1	95.0	d__Archaea;p__Methanobacteriota;c__Methanobacteria;o__Methanobacteriales;f__Methanobacteriaceae;g__Methanobrevibacter;s__Methanobrevibacter ruminantium	100.0	1.0	N/A	N/A	N/A	N/A	N/A	N/A	ani_screen	classification based on ANI only	GCA_900321995.1, s__Methanobrevibacter sp900321995, 95.0, 80.9, 0.7; GCF_900114585.1, s__Methanobrevibacter olleyae, 95.0, 79.96, 0.55	N/A	N/A	N/A	N/A
    genome_2	d__Archaea;p__Thermoplasmatota;c__Thermoplasmata;o__Methanomassiliicoccales;f__Methanomethylophilaceae;g__VadinCA11;s__VadinCA11 sp002498365	GCA_002498365.1	95.0	d__Archaea;p__Thermoplasmatota;c__Thermoplasmata;o__Methanomassiliicoccales;f__Methanomethylophilaceae;g__VadinCA11;s__VadinCA11 sp002498365	99.16	0.94	GCA_002498365.1	95.0	d__Archaea;p__Thermoplasmatota;c__Thermoplasmata;o__Methanomassiliicoccales;f__Methanomethylophilaceae;g__VadinCA11;s__VadinCA11 sp002498365	99.16	0.94	d__Archaea;p__Thermoplasmatota;c__Thermoplasmata;o__Methanomassiliicoccales;f__Methanomethylophilaceae;g__VadinCA11;s__	taxonomic classification defined by topology and ANI	topological placement and ANI have congruent species assignments	GCA_002505345.1, s__VadinCA11 sp002505345, 95.0, 89.92, 0.89; GCA_002509405.1, s__VadinCA11 sp002509405, 95.0, 88.13, 0.89	87.1	11	N/A	N/A
    genome_3	d__Archaea;p__Thermoproteota;c__Nitrososphaeria;o__Nitrososphaerales;f__Nitrosopumilaceae;g__;s__	N/A	N/A	N/A	N/A	N/A	N/A	N/A	N/A	N/A	N/A	d__Archaea;p__Thermoproteota;c__Nitrososphaeria;o__Nitrososphaerales;f__Nitrosopumilaceae;g__;s__	taxonomic novelty determined using RED	N/A	N/A	85.47	11	0.89523	N/A
    genome_4	Unclassified	N/A	N/A	N/A	N/A	N/A	N/A	N/A	N/A	N/A	N/A	N/A	N/A	N/A	N/A	N/A	N/A	N/A	No bacterial or archaeal marker
