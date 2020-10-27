.. _files/summary.tsv:


summary.tsv
===========

Classifications provided by the GTDB-Tk are in the files \<prefix>.bac120.summary.tsv and \<prefix>.ar122.summary.tsv for bacterial and archaeal genomes, respectively. These are tab separated files with the following columns:

* user_genome: Unique identifier of query genome taken from the FASTA file of the genome.
* classification: GTDB taxonomy string inferred by the GTDB-Tk. An unassigned species (i.e., ``s__``) indicates that the query genome is either i) placed outside a named genus or ii) the ANI to the closest intra-genus reference genome with an AF >=0.65 is not within the species-specific ANI circumscription radius.
* fastani_reference: indicates the accession number of the reference genome (species) to which a user genome was assigned based on ANI and AF. ANI values are only calculated when a query genome is placed within a defined genus and are evaluated for all reference genomes in that genus.
* fastani_reference_radius: indicates the species-specific ANI circumscription radius of the reference genomes used to determine if a query genome should be classified to the same species as the reference.
* fastani_taxonomy: indicates the GTDB taxonomy of the above reference genome.
* fastani_ani: indicates the ANI between the query and above reference genome.
* fastani_af: indicates the AF between the query and above reference genome.
* closest_placement_reference: indicates the accession number of the reference genome when a genome is placed on a terminal branch. 
* closest_placement_taxonomy: indicates the GTDB taxonomy of the above reference genome.
* closest_placement_ani: indicates the ANI between the query and above reference genome.
* closest_placement_af: indicates the AF between the query and above reference genome.
* pplacer_taxonomy: indicates the pplacer taxonomy of the query genome.
* classification_method:	indicates the rule used to classify the genome. This field will be one of: i) ANI, indicating a species assignement was based solely on the calculated ANI and AF with a reference genome; ii) ANI/Placement, indicating a species assignment was made based on both ANI and the placement of the genome in the reference tree; iii) taxonomic classification fully defined by topology, indicating that the classification could be determine based solely on the genome's position in the reference tree; or iv) taxonomic novelty determined using RED, indicating that the relative evolutionary divergence (RED) and placement of the genome in the reference tree were used to determine the classification.
* note: provides additional information regarding the classification of the genome. Currently this field is only filled out when a species determination is made and indicates if the placement of the genome in the reference tree and closest reference according to ANI/AF are the same (congruent) or different (incongruent). 
* other_related_references: lists up to the 100 closest reference genomes based on ANI. ANI calculations are only performed between a query genome and reference genomes in the same genus.
* msa_percent: indicates the percentage of the MSA spanned by the genome (i.e. percentage of columns with an amino acid).
* red_value: indicates, when required, the relative evolutionary divergence (RED) for a query genome. RED is not calculated when a query genome can be classified based on ANI.
* warnings: indicates unusual characteristics of the query genome that may impact the taxonomic assignment.

Produced by
-----------

 * :ref:`commands/classify`
 * :ref:`commands/classify_wf`

Example
-------

.. code-block:: text

    user_genome	classification	fastani_reference	fastani_reference_radius	fastani_taxonomy	fastani_ani	fastani_af	closest_placement_reference	closest_placement_taxonomy	closest_placement_ani	closest_placement_af	pplacer_taxonomy	classification_method	note	other_related_references(genome_id,species_name,radius,ANI,AF)	aa_percent	translation_table	red_value	warnings
    genome_2	d__Archaea;p__Thermoplasmatota;c__Thermoplasmata;o__Methanomassiliicoccales;f__Methanomethylophilaceae;g__VadinCA11;s__VadinCA11 sp002498365	GCA_002498365.1	95.0	d__Archaea;p__Thermoplasmatota;c__Thermoplasmata;o__Methanomassiliicoccales;f__Methanomethylophilaceae;g__VadinCA11;s__VadinCA11 sp002498365	99.16	0.94	GCA_002498365.1	d__Archaea;p__Thermoplasmatota;c__Thermoplasmata;o__Methanomassiliicoccales;f__Methanomethylophilaceae;g__VadinCA11;s__VadinCA11 sp002498365	99.16	0.94	d__Archaea;p__Thermoplasmatota;c__Thermoplasmata;o__Methanomassiliicoccales;f__Methanomethylophilaceae;g__VadinCA11;s__	ANI/Placement	topological placement and ANI have congruent species assignments	GCA_002505345.1, s__VadinCA11 sp002505345, 95.0, 89.92, 0.89; GCA_002509405.1, s__VadinCA11 sp002509405, 95.0, 88.13, 0.89	87.1	11	N/A	N/A
    genome_3	d__Archaea;p__Thermoplasmatota;c__Thermoplasmata;o__Methanomassiliicoccales;f__Methanomethylophilaceae;g__VadinCA11;s__VadinCA11 sp002498365	GCA_002498365.1	95.0	d__Archaea;p__Thermoplasmatota;c__Thermoplasmata;o__Methanomassiliicoccales;f__Methanomethylophilaceae;g__VadinCA11;s__VadinCA11 sp002498365	95.33	0.87	GCA_002498365.1	d__Archaea;p__Thermoplasmatota;c__Thermoplasmata;o__Methanomassiliicoccales;f__Methanomethylophilaceae;g__VadinCA11;s__VadinCA11 sp002498365	95.33	0.87	d__Archaea;p__Thermoplasmatota;c__Thermoplasmata;o__Methanomassiliicoccales;f__Methanomethylophilaceae;g__VadinCA11;s__	ANI/Placement	topological placement and ANI have congruent species assignments	GCA_002505345.1, s__VadinCA11 sp002505345, 95.0, 94.26, 0.87; GCA_002509405.1, s__VadinCA11 sp002509405, 95.0, 90.74, 0.77	73.07	11	N/A	N/A
    genome_1	d__Archaea;p__Euryarchaeota;c__Methanobacteria;o__Methanobacteriales;f__Methanobacteriaceae;g__Methanobrevibacter;s__Methanobrevibacter ruminantium	GCF_000024185.1	95.0	d__Archaea;p__Euryarchaeota;c__Methanobacteria;o__Methanobacteriales;f__Methanobacteriaceae;g__Methanobrevibacter;s__Methanobrevibacter ruminantium	100.0	1.0	GCF_000024185.1	d__Archaea;p__Euryarchaeota;c__Methanobacteria;o__Methanobacteriales;f__Methanobacteriaceae;g__Methanobrevibacter;s__Methanobrevibacter ruminantium	100.0	1.0	d__Archaea;p__Euryarchaeota;c__Methanobacteria;o__Methanobacteriales;f__Methanobacteriaceae;g__Methanobrevibacter;s__	ANI/Placement	topological placement and ANI have congruent species assignments	GCA_900321995.1, s__Methanobrevibacter sp900321995, 95.0, 80.9, 0.7; GCF_900114585.1, s__Methanobrevibacter olleyae, 95.0, 79.96, 0.55; GCA_900314635.1, s__Methanobrevibacter sp900314635, 95.0, 78.45, 0.3	97.09	11	N/A	N/A


