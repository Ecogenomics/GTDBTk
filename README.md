# GTDB-Tk

[![version status](https://img.shields.io/pypi/v/gtdbtk.svg)](https://pypi.python.org/pypi/gtdbtk)
[![Bioconda](https://img.shields.io/conda/vn/bioconda/gtdbtk.svg?color=green)](https://anaconda.org/bioconda/gtdbtk)
[![Downloads](https://pepy.tech/badge/gtdbtk/month)](https://pepy.tech/project/gtdbtk)

GTDB-Tk is a software toolkit for assigning objective taxonomic classifications to bacterial and archaeal genomes. It is designed to work with recent advances that allow hundreds or thousands of metagenome-assembled genomes (MAGs) to be obtained directly from environmental samples. It can also be applied to isolate and single-cell genomes. The GTDB-Tk is open source and released under the [GNU General Public License (Version 3)](https://www.gnu.org/licenses/gpl-3.0.en.html).

Notifications about GTDB-Tk releases will be available through the GTDB Twitter account (https://twitter.com/ace_gtdb).

Please visit the GTDB-Tk [github page](https://github.com/Ecogenomics/GTDBTk) for the latest updates about the software.


* [Announcements](#announcements)
* [Hardware requirements](#hardware-requirements)
* [Dependencies](#dependencies)
  * [Python libraries](#python-libraries)
  * [Third-party software](#third-party-software)
  * [GTDB-Tk reference data](#gtdb-tk-reference-data)
* [Installation](#installation)
  * [pip installation](#pip-installation)
  * [Bioconda installation](#bioconda-installation)
  * [Testing installation](#testing-installation)
* [FAQ](docs/faq.md)
* [Quick start](#quick-start)
* [Classify workflow](#classify-workflow)
* [Validating species assignments with ANI](#validating-species-assignments-with-average-nucleotide-identity)
* [Classification summary file](#classification-summary-file)
* [<i>De novo</i> workflow](#de-novo-workflow)
* [Individual steps](#individual-steps)
* [References](#references)

## Announcements

**Note (July 12, 2019)**:
* GTDB-Tk v0.3.2 has been released (**we recommend all users update to this version**)
    * FastANI calculations are more robust.
    * Optimisation of RED calculations.
    * Improved output messages when errors are encountered.

**Note (July 08, 2019)**:
* GTDB-Tk v0.3.1 has been released:
  * Pplacer taxonomy is now available in the summary file.
  * FastANI species assignment will be selected over phylogenetic placement (Topology case).
  
**Note (June 21, 2019)**:
* GTDB-Tk v0.3.0 has been released:
  * Best translation table displayed in summary file.
  * GTDB-Tk now supports gzipped genomes as inputs (--extension .gz).
  * By default, GTDB-Tk uses precalculated RED values.
  * New option to recalculate RED value during classify step (--recalculate_red).
  * New option to export the untrimmed reference MSA files.
  * New option to skip_trimming during align step.
  * New option to use a custom taxonomy file when rooting a tree.
  * New [FAQ](docs/faq.md) page available.
  * New output structure.
  * This version requires a new version of the GTDB-Tk data package (gtdbtk_r89_data.tar.gz) available [here](https://data.ace.uq.edu.au/public/gtdb/data/releases/release89/89.0/)
  
[Previous announcements](docs/announcements.md)


## Hardware requirements
- ~100Gb of memory to run
- ~27Gb of storage
- ~1 hour per 1,000 genomes when using 64 CPUs

## Dependencies

### Python libraries

GTDB-Tk is designed for Python 2.7 and requires the following Python libraries:
* [dendropy](http://dendropy.org/)  >=4.1.0: Python library for phylogenetics.
* [future](https://python-future.org/index.html) >= 0.15.2: Clean single-source support for Python 3 and 2
* [SciPy Stack](https://www.scipy.org/install.html): at least the Matplotlib, NumPy, and SciPy libraries.

Dendropy and future will be installed as part of GTDB-Tk when installing via pip (see below). The **SciPy Stack** must be installed separately.

### Third-party software

GTDB-Tk makes use of the following 3rd party dependencies and assumes these are on your system path:
* [Prodigal](http://compbio.ornl.gov/prodigal/) >= 2.6.2: Hyatt D, et al. 2012. Gene and translation initiation site prediction in metagenomic sequences. <i>Bioinformatics</i>, 28, 2223-2230.
* [HMMER](http://hmmer.org/) >= 3.1: Eddy SR. 2011. Accelerated profile HMM searches. <i>PLoS Comp. Biol.</i>, 7, e1002195.
* [pplacer](http://matsen.fhcrc.org/pplacer/) >= 1.1: Matsen F, et al. 2010. pplacer: linear time maximum-likelihood and Bayesian phylogenetic placement of sequences onto a fixed reference tree. <i>BMC Bioinformatics</i>, 11, 538.
* [FastANI](https://github.com/ParBLiSS/FastANI) >= 1.0: Jain C, et al. 2018. High-throughput ANI analysis of 90K prokaryotic genomes reveals clear species boundaries. <i>Nature Communication</i>, 5114.
* [FastTree](http://www.microbesonline.org/fasttree/) >= 2.1.9: Price MN, et al. 2010 FastTree 2 -- Approximately Maximum-Likelihood Trees for Large Alignments. <i>PLoS ONE</i>, 5, e9490.

Please cite these tools if you use GTDB-Tk in your work.

### GTDB-Tk reference data

GTDB-Tk requires ~27G+ of external data that need to be downloaded and unarchived:
```
wget https://data.ace.uq.edu.au/public/gtdb/data/releases/release89/89.0/gtdbtk_r89_data.tar.gz
tar xvzf gtdbtk_r89_data.tar.gz
```

Reference data for prior releases of GTDB-Tk are available at:
```
wget https://data.ace.uq.edu.au/public/gtdbtk
```

## Installation

### pip installation

Once dependencies are installed, GTDB-Tk can be installed using [pip](https://pypi.python.org/pypi/gtdbtk):
```
> pip install gtdbtk
```
GTDB-Tk requires an environmental variable named GTDBTK_DATA_PATH to be set to the directory containing the data downloaded from https://data.ace.uq.edu.au/public/gtdbtk/.
```
export GTDBTK_DATA_PATH=/path/to/release/package/
```
Alternatively, you can permanently add this variable to your .bash_profile as described [here](https://unix.stackexchange.com/questions/26047/how-to-correctly-add-a-path-to-path).

You may also wish to add the GTDB-Tk binary file to your .bash_profile:
```
echo 'alias gtdbtk="python ~/.local/bin/gtdbtk"' >> ~/.bashrc
```

### Bioconda installation

A Bioconda recipe is available as an alternative installation source (thanks to [Natasha](https://github.com/npavlovikj) for starting it). You can find the recipe here: https://anaconda.org/bioconda/gtdbtk

The GTDB-Tk reference data is not available as a part of the recipe, but it can be automatically downloaded to the conda package by executing the command: `download-db.sh`


### Testing installation

You can test your GTDB-Tk installation by running:
```
gtdbtk test --out_dir <path>
```
This applies the classify workflow (classify_wf) to three archaeal genomes. It creates two folders in the output directory: the result directory and the genomes directory which is used as input to GTDB-Tk.

## Quick start

The functionality provided by GTDB-Tk can be accessed through the help menu:
```
> gtdbtk -h
```

Usage information about each method can also be accessed through their specific help menu, e.g.:
```
> gtdbtk classify_wf -h
```

## Classify workflow

The classify workflow consists of three steps: *identify*, *align*, and *classify*. The *identify* step calls genes using [Prodigal](http://compbio.ornl.gov/prodigal/), and uses HMM models and the [HMMER](http://hmmer.org/) package to identify the 120 bacterial and 122 archaeal marker genes used for phylogenetic inference. Multiple sequence alignments (MSA) are obtained by aligning marker genes to their respective HMM model. The *align* step concatenates the aligned marker genes and filters the concatenated MSA to approximately 5,000 amino acids. Finally, the *classify* step uses [pplacer](http://matsen.fhcrc.org/pplacer/) to find the maximum-likelihood placement of each genome in the GTDB-Tk reference tree. GTDB-Tk classifies each genome based on its placement in the reference tree, its relative evolutionary divergence, and/or average nucleotide identity (ANI) to reference genomes.
 
The classify workflow can be run as follows:
```
> gtdbtk classify_wf --genome_dir <my_genomes> --out_dir <output_dir>
```
This will process all genomes in <my_genomes> using both bacterial and archaeal marker sets and place the results in <output_dir>. Genomes must be in FASTA format. The location of genomes can also be specified using a batch file with the --batchfile flag. The batch file is a two column file indicating the location of each genome and the desired genome identifier (i.e., a Newick compatible alphanumeric string). These fields must be separated by a tab.

The workflow supports several optional flags, including:
* min_perc_aa: allows filtering of genomes below a specified percentage of amino acids in the MSA
* cpus: maximum number of CPUs to use

For other flags please consult the command line interface.

Here is an example run of this workflow:
```
> gtdbtk classify_wf --cpus 24 --genome_dir ./my_genomes --out_dir gtdbtk_output
```

The taxonomic classification of each bacterial and archaeal genome is contained in the \<prefix\>.bac120.summary.tsv and \<prefix\>.ar122.summary.tsv output files.

##### Additional output files 

Each step of the classify workflow generates a number of files that can be consulted for additional information about the processed genomes.

Identify step:
* identify/\<prefix\>.bac120.markers_summary.tsv: summary of unique, duplicated, and missing markers within the 120 bacterial marker set for each submitted genome.
* identify/\<prefix\>.ar122.markers_summary.tsv: analogous to the above file, but for the 122 archaeal marker set.
* identify/\<prefix\>.translation_table_summary.tsv: The predicted [translation table](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi) used for gene calling for each genome.
* identify/intermediate_results/marker_genes/: contains individual genome results for gene calling using Prodigal and gene identification based on TIGRFAM and Pfam HMMs.

Align step:
* align/\<prefix\>.\[bac120/ar122\].user_msa.fasta: FASTA file containing MSA of the submitted genomes.
* align/\<prefix\>.\[bac120/ar122\].msa.fasta: FASTA file containing MSA of submitted and reference genomes.
* align/\<prefix\>.\[bac120/ar122\].filtered.tsv: list of genomes with an insufficient number of amino acids in MSA.
* align/intermediate_results/\<prefix\>.\[bac120/ar122\].marker_info.tsv: Markers used in generation of the concatenated MSA and the order in which they were applied.

Classify step:
* classify/\<prefix>.\[bac120/ar122\].classify.tree: reference tree in Newick format containing query genomes placed with pplacer.
* classify/\<prefix>.\[bac120/ar122\].summary.tsv: classification of query genomes based on their placement in the reference tree, relative evolutionary divergence, and ANI to reference genomes. This is the primary output of the GTDB-Tk and contains the taxonomic classification we recommend plus additional information regarding the criteria used to assign taxonomy (see below).
* classify/intermediate_results/\<prefix>.\[bac120/ar122\].classification_pplacer.tsv: classification of query genomes based only on their placement in the reference tree.
* classify/intermediate_results/\<prefix>.\[bac120/ar122\].red_dictionary.tsv: median RED values for taxonomic ranks.
* classify/intermediate_results/pplacer/: Output information generated by pplacer. 

## Validating species assignments with average nucleotide identity

The GTDB-Tk uses [FastANI](https://github.com/ParBLiSS/FastANI) to estimate the ANI between genomes. A query genome is classified as belonging to the same species as a reference genome if the ANI between the genomes is within the species ANI circumscription radius (typically, 95%) and the alignment fraction (AF) is >= 0.65 . In exceptional circumstances, the phylogenetic placement of a query genome may not support the species assignment. GTDB r89 will strictly use ANI to circumscribe species and GTDB-Tk follows this methodology.

## Classification summary file 

Classifications provided by the GTDB-Tk are in the files \<prefix>.bac120.summary.tsv and \<prefix>.ar122.summary.tsv for bacterial and archaeal genomes, respectively. These are tab separated files with the following columns:

* user_genome: Unique identifier of query genome taken from the FASTA file of the genome.
* classification: GTDB taxonomy string inferred by the GTDB-Tk. An unassigned species (i.e., s__) indicates that the query genome is either i) placed outside a named genus or ii) has an ANI <95% to all reference genomes within the same genus as the query genome.
* fastani_reference: Indicates the accession number of the closest reference genome as determine by ANI and AF. This genome is used along with the placement of the genome in the reference tree to determine the species assignment on the genome. ANI values are only calculated when a query genome is placed within a defined genus and are calculated for all reference genomes in that genus.
* fastani_reference_radius: Indicates the ANI threshold of the reference genomes used to determine if a query genome should be classified to the same species as the reference.
* fastani_taxonomy: Indicates the GTDB taxonomy of the above reference genome.
* fastani_ani: Indicates the ANI between the query and above reference genome.
* fastani_af: Indicates the alignment fraction (AF) between the query and above reference genome.
* closest_placement_reference: Indicates the accession number of the reference genome when a genome is placed on a terminal branch. This genome is used along with the ANI information to determine the species assignment on the genome.
* closest_placement_taxonomy: Indicates the GTDB taxonomy of the above reference genome.
* closest_placement_ani: Indicates the ANI between the query and above reference genome.
* closest_placement_af: Indicates the AF between the query and above reference genome.
* pplacer_taxonomy: Indicates the pplacer taxonomy of the query genome.
* classification_method:	Indicates the rule used to classify the genome. This field will be one of: i) ANI, indicating a species assignement was based solely on the calculated ANI with a reference genome; ii) ANI/Placement, indicating a species assignment was made based on both ANI and the placement of the genome in the reference tree; iii) taxonomic classification fully defined by topology, indicating that the classification could be determine based solely on the genome's position in the reference tree; or iv) taxonomic novelty determined using RED, indicating that the relative evolutionary divergence (RED) and placement of the genome in the reference tree were used to determine the classification.
* note: Provides additional information regarding the classification of the genome. Currently this field is only filled out when a species determination must be made and indicates that the placement of the genome and closest reference according to ANI are either the same (congruent) or different (incongruent). 
* other_related_references: Lists up to the top 100 closest reference genomes based on ANI. ANI calculations are only performed between a query genome and reference genomes in the same genus.
* aa_percent: Indicates the percentage of the MSA spanned by the genome (i.e. percentage of columns with an amino acid).
* red_value: Indicates, when required, the relative evolutionary divergence (RED) for a query genome. RED is not calculated when a query genome can be classified based on ANI.
* warnings: Indicates unusual characteristics of the query genome that may impact the taxonomic assignment

## De novo workflow

**under active development - decorate step not yet implemented**

The *de novo* workflow infers new bacterial and archaeal trees containing all user supplied and GTDB-Tk reference genomes. The classify workflow is recommended for obtaining taxonomic classifications, and this workflow only recommended if a *de novo* domain-specific trees are desired. This workflow consists of five steps: *identify*, *align*, *infer*, *root*, and *decorate* (not yet implemented). The *identify* and *align* steps are the same as in the classify workflow. The *infer* step uses [FastTree](http://www.microbesonline.org/fasttree/) with the WAG+GAMMA models to calculate independent, *de novo* bacterial and archaeal trees. These trees can then be rooted using a user specified outgroup and decorated with the GTDB taxonomy. 

The *de novo* workflow can be run as follows:
```
> gtdbtk de_novo_wf --genome_dir <my_genomes> --<marker_set> --outgroup_taxon <outgroup> --out_dir <output_dir>
```
This will process all genomes in <my_genomes> using the specified marker set and place the results in <output_dir>. Only genomes previously identified as being bacterial (archaeal) should be included when using the bacterial (archaeal) marker set. The tree will be rooted with the <outgroup> taxon (typically a phylum in the domain-specific tree) as required for correct decoration of the tree. In general, we suggest the resulting tree be treated as unrooted when interpreting results. Identical to the classify workflow, the location of genomes can also be specified using a batch file with the --batchfile flag.

The workflow supports several optional flags, including:
* cpus: maximum number of CPUs to use
* min_perc_aa: filter genomes with an insufficient percentage of AA in the MSA (default: 50)
* taxa_filter: filter genomes to taxa within specific taxonomic groups
* prot_model: protein substitution model for tree inference (LG or WAG; default: WAG)

For other flags please consult the command line interface.

Here is an example run of this workflow:
```
> gtdbtk de_novo_wf --genome_dir ./genomes --bac120_ms --outgroup_taxon p__Chloroflexota --taxa_filter p__Firmicutes --out_dir de_novo_output
```

## Individual steps

All steps comprising the classify and <i>de novo</i> workflows can be run independently if desired. Please consult the command line interface for specific details on running each of these steps.

## References

A manuscript describing the GTDB-Tk is currently being prepared:

* Chaumeil PA, Mussig AJ, Hugenholtz P, Parks DH. 2019. GTDB-Tk: A toolkit to classify genomes with the Genome Taxonomy Database. \<in prep\>.

In the meantime, if you find the GTDB-Tk useful please cite this GitHub page or the GTDB taxonomy:

* Parks DH, et al. 2018. [A standardized bacterial taxonomy based on genome phylogeny substantially revises the tree of life](https://www.nature.com/articles/nbt.4229). <i>Nat. Biotechnol.</i>, http://dx.doi.org/10.1038/nbt.4229
 
 We also strongly encourage you to cite the following 3rd party dependencies:

* Matsen FA, Kodner RB, Armbrust EV. 2010. [pplacer: linear time maximum-likelihood and Bayesian phylogenetic placement of sequences onto a fixed reference tree](https://www.ncbi.nlm.nih.gov/pubmed/21034504). <i>BMC Bioinformatics</i>, 11:538.
* Jain C, et al. 2019. [High-throughput ANI Analysis of 90K Prokaryotic Genomes Reveals Clear Species Boundaries](https://www.nature.com/articles/s41467-018-07641-9). <i>Nat. Communications</i>, doi: 10.1038/s41467-018-07641-9.
* Hyatt D, et al. 2010. [Prodigal: prokaryotic gene recognition and translation initiation site identification](https://www.ncbi.nlm.nih.gov/pubmed/20211023). <i>BMC Bioinformatics</i>, 11:119. doi: 10.1186/1471-2105-11-119.
* Price MN, Dehal PS, Arkin AP. [FastTree 2 - Approximately Maximum-Likelihood Trees for Large Alignments](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2835736/). <i>PLoS One</i>, 5, e9490.
* Eddy SR. 2011. [Accelerated profile HMM searches](https://www.ncbi.nlm.nih.gov/pubmed/22039361). <i>PLOS Comp. Biol.</i>, 7:e1002195.

## Copyright

Copyright 2017 Pierre-Alain Chaumeil. See LICENSE for further details.
