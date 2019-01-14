# GTDB-Tk

[![version status](https://img.shields.io/pypi/v/gtdbtk.svg)](https://pypi.python.org/pypi/gtdbtk)

GTDB-Tk is a software toolkit for assigning objective taxonomic classifications to bacterial and archaeal genomes. It is computationally 
efficient and designed to work with recent advances that allow hundreds or thousands of metagenome-assembled genomes (MAGs) to be obtained directly from environmental samples. It can also be applied to isolate and single-cell genomes. The GTDB-Tk is open source and released under the GNU General Public License (Version 3).

GTDB-Tk is **under active development and validation**. Please independently confirm the GTDB-Tk predictions by manually inspecting the tree and bring any discrepencies to our attention. Notifications about GTDB-Tk releases will be available through the GTDB Twitter account (https://twitter.com/ace_gtdb).


* [Announcements](#announcements)
* [Hardware requirements](#hardware-requirements)
* [Dependencies](#dependencies)
  * [Python libraries](#python-libraries)
  * [Third-party software](#third-party-software)
  * [Perl modules](#perl-modules)
  * [GTDB-Tk reference data](#gtdb-tk-reference-data)
* [pip installation](#pip-installation)
* [Bioconda installation](#bioconda-installation)
* [Quick start](#quick-start)
* [Classify workflow](#classify-workflow)
* [Validating species assignments with ANI](#validating-species-assignments-with-average-nucleotide-identity)
* [Classification summary file](#classification-summary-file)
* [<i>De novo</i> workflow](#de-novo-workflow)
* [Individual steps](#individual-steps)
* [References](#references)

## Announcements
**Note (Sept. 20, 2018)**:
- GTDB-Tk v0.1.3 has been released and addresses an issue with species assignments based on the placement of genomes in the reference tree. This impacted species assignment when submitting multiple closely related genomes. Species assignments reported by ANI were not impacted.
- We recommend all users update to this version.

**Note (Aug. 30, 2018)**:
- A new version of the data (release 86) is available under [this link](https://data.ace.uq.edu.au/public/gtdbtk/release_86/).
- This new version is required to run GTDB-Tk v0.1.0+

## Hardware requirements
- ~90Gb of memory to run
- ~25Gb of storage
- ~1 hour per 1,000 genomes when using 64 CPUs

## Dependencies

### Python libraries

GTDB-Tk requires the following Python libraries:
* [jinja2](http://jinja.pocoo.org/) >=2.7.3: a full featured template engine for Python.
* [mpld3](http://mpld3.github.io/) >=0.2: D3 viewer for Matplotlib.
* [biolib](https://github.com/dparks1134/biolib) >=0.0.44: Python package for common tasks in bioinformatic.
* [dendropy](http://dendropy.org/)  >=4.1.0: Python library for phylogenetics.
* [SciPy Stack](https://www.scipy.org/install.html): at least the Matplotlib, NumPy, and SciPy libraries.

Jinja2, mpld3, dendropy and biolib will be installed as part of GTDB-Tk when installing it via pip (see below). The **SciPy Stack** must be installed separately.

### Third-party software

GTDB-Tk makes use of the following 3rd party dependencies and assumes these are on your system path:
* [Prodigal](http://prodigal.ornl.gov/) >= 2.6.2: Hyatt D, et al. 2012. Gene and translation initiation site prediction in metagenomic sequences. <i>Bioinformatics</i>, 28, 2223-2230.
* [HMMER](http://http://hmmer.org/) >= 3.1: Eddy SR. 2011. Accelerated profile HMM searches. <i>PLoS Comp. Biol.</i>, 7, e1002195.
* [pplacer](http://matsen.fhcrc.org/pplacer/) >= 1.1: Matsen F, et al. 2010. pplacer: linear time maximum-likelihood and Bayesian phylogenetic placement of sequences onto a fixed reference tree. <i>BMC Bioinformatics</i>, 11, 538.
* [FastANI](https://github.com/ParBLiSS/
) >= 1.0: Jain C, et al. 2017. High-throughput ANI Analysis of 90K Prokaryotic Genomes Reveals Clear Species Boundaries.<i>bioRxiv.</i> 256800.
* [FastTree](http://www.microbesonline.org/fasttree/) >= 2.1.9: Price MN, et al. 2010 FastTree 2 -- Approximately Maximum-Likelihood Trees for Large Alignments. <i>PLoS ONE</i>, 5, e9490.

### Perl modules

GTDB-Tk assumes the Python 2.7.x and Perl interpreters are on your system path.    
<u>**note:**</u>  Perl interpreter requires Moose, Bundle::BioPerl and IPC::Run modules. You can install those modules using CPAN:
```
perl -MCPAN -e"install Moose"
perl -MCPAN -e"install IPC::Run"
perl -MCPAN -e"install Bundle::BioPerl"
```
If ```perl -MCPAN -e"install Bundle::BioPerl"``` does not run on your server, please install BioPerl using the steps under [this link](https://bioperl.org/INSTALL.html).

Make sure that the installed Perl modules (.pm) paths are part of the @inc variable. If not, the PERL5LIB (or PERLIB) environment variable need to be updated the same way the PATH environment variable is updated. Every directory listed in this variable will be added to @inc, e.g.:
```
export PERL5LIB="$PERL5LIB:/path/to/moose/module:/path/to/ipc/module:/path/to/bioperl/module"
```

### GTDB-Tk reference data

GTDB-Tk requires ~25G+ of external data that need to be downloaded and unarchived:
```
wget https://data.ace.uq.edu.au/public/gtdbtk/release_86/gtdbtk_r86_data.tar.gz
tar xvzf gtdbtk_r86_archived_data.tar.gz
```

Reference data for prior releases of GTDB-Tk are available at:
```
wget https://data.ace.uq.edu.au/public/gtdbtk
```

## pip installation

Once dependencies are installed, GTDB-Tk can be installed using [pip](https://pypi.python.org/pypi/gtdbtk):
```
> pip install gtdbtk
```

GTDB-Tk requires a config file. In the Python lib/site-packages directory, go to the gtdbtk directory and setup this config file:
```
cd config
cp config_template.py config.py
```
Edit the config.py file and modify the GENERIC_PATH variables so it point to the directory containing the data downloaded from https://data.ace.uq.edu.au/public/gtdbtk/. Make sure the variable finishes with a slash '/'.

## Bioconda installation

A Bioconda receipe has been put together by [Natasha](https://github.com/npavlovikj) (thanks!). You can find the receipe at:
https://anaconda.org/bioconda/gtdbtk

The download of the GTDB-Tk reference data is not part of the recipe, but there is a "download-db.sh" script that does that when run from the conda environment.

## Quick start

The functionality provided by GTDB-Tk can be accessed through the help menu:
```
> gtdbtk -h
```

Usage information about each methods can also be accessed through their species help menu, e.g.:
```
> gtdbtk classify_wf -h
```

## Classify workflow

The classify workflow consists of three steps: *identify*, *align*, and *classify*. The *identify* step calls genes using [Prodigal](http://prodigal.ornl.gov/), and HMM models and the [HMMER](http://http://hmmer.org/) package to identify the 120 bacterial and 122 archaeal marker genes used for phylogenetic inference. Consistent alignments are obtained by aligning marker genes to their respective HMM model. The *align* step concatenates the aligned marker genes and filters the concatenated multiple sequence alignment (MSA) to approximately 5,000 amino acids. Finally, the *classify* step uses [pplacer](http://matsen.fhcrc.org/pplacer/) to find the maximum-likelihood placement of each genome in the GTDB-Tk reference tree. GTDB-Tk classifies each genome based on its placement in the reference tree, its relative evolutionary divergence, and ANI to reference genomes.
 
The classify workflow can be run as follows:
```
> gtdbtk classify_wf --genome_dir <my_genomes> --out_dir <output_dir>
```
This will process all genomes in <my_genomes> using both bacterial and archaeal marker sets and place the results in <output_dir>. Genomes must be in FASTA format. The location of genomes can also be specified using a batch file with the --batchfile flag. The batch file is a two column file indicating the location of each genome and the desired genome identifier (i.e., a Newick compatible alphanumeric string). These fields must be seperated by a tab.

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
* \<prefix\>_bac120_markers_summary.tsv: summary of unique, duplicated, and missing markers within the 120 bacterial marker set for each submitted genome
* \<prefix\>_ar122_markers_summary.tsv: analogous to the above file, but for the 122 archaeal marker set
* marker_genes directory: contains individual genome results for gene calling using Prodigal and gene identification based on TIGRFAM and Pfam HMMs

Align step:
* \<prefix\>.user_msa.fasta: FASTA file containing MSA of the submitted genomes
* \<prefix\>.msa.fasta: FASTA file containing MSA of submitted and reference genomes
* \<prefix\>.filtered.tsv: list of genomes with an insufficient number of amino acids in MSA

Classify step:
 * \<prefix>.summary.tsv: classification of user genomes based on their placement in the reference tree, relative evolutionary divergence, and ANI to reference genomes. This is the primary output of the GTDB-Tk and contains the taxonomic classification we recommend plus additional information regarding the criteria used to assign genomes (see below)
* \<prefix>.classification_pplacer.tsv: classification of user genomes based only on their pplacer placement in the reference tree
* \<prefix>.classify.tree: reference tree in Newick format containing user genomes placed with pplacer
* \<prefix>.red_dictionary.tsv: median RED values for taxonomic ranks

## Validating species assignments with average nucleotide identity

The GTDB-Tk uses [FastANI](https://github.com/ParBLiSS/FastANI) to estimate the average nucleotide identity (ANI) between genomes. Reference genomes with an ANI >= 95% to a query genome are used to validate species assignments.

## Classification summary file 

Classifications provided by the GTDB-Tk are in the files \<prefix>.bac120.summary.tsv and \<prefix>.ar122.summary.tsv for bacterial and archaeal genomes, respectively. These are tab separated files with the following columns:

* user_genome: Unique identifer of genome taken from the FASTA file for the genome.
* classification: GTDB taxonomy string inferred by the GTDB-Tk. The GTDB does not provide species assignments for all reference genomes (e.g., taxonomic groups composed entirely of metagenome-assembled genomes). In such cases, species assignments will reflect the accession number of the reference genome (e.g., s__GCA_001940855.1). A unassigned species (i.e., s__) indicates that either the genome represents a novel species or that a species assignment could not be reliably established as indicated by the following fields.
* fastani_reference: Indicates the accession number of the closest reference genome as determine by ANI. This genome is used along with the placement of the genome in the reference tree to determine the species assignment on the genome. ANI values are only calculated when a query genome is placed within a defined genus and are calculated for all reference genomes in the genus.
* fastani_taxonomy: Indicates the GTDB taxonomy of the above reference genome.
* fastani_ani: Indicates the ANI between the query and above reference genome.
* fastani_af: Indicates the AF between the query and above reference genome.
* closest_placement_reference: Indicates the accession number of the reference genome when a genome is placed on a terminal branch. This genome is used along with the ANI information to determine the species assignment on the genome.
* closest_placement_taxonomy: Indicates the GTDB taxonomy of the above reference genome.
* closest_placement_ani: Indicates the ANI between the query and above reference genome.
* closest_placement_af: Indicates the AF between the query and above reference genome.
* classification_method:	Indicates the rule used to classify the genome. This field will be one of: i) ANI/Placement, indicating a species assignment was made based on both the calculate ANI and placement of the genome in the reference tree; ii) taxonomic classification fully defined by topology, indicating that the classification could be determine based solely on the genome's position in the reference tree; or iii) taxonomic novelty determined using RED, indicating that the relative evolutionary divergence (RED) and placement of the genome in the reference tree were used to determine the classification.
* note: Provides additional information regarding the classification of the genome. Currently this field is only filled out when a species determination must be made and indicates that the placement of the genome and closest reference according to ANI are either the same (congruent) or different (incongruent). The GTDB-Tk will leave the species field empty (i.e., s__) when the two methods are incongruent.
* other_related_references: Lists the top 100 closest reference genomes.

## De novo workflow

**under active development**

The *de novo* workflow infers a new tree containing all user supplied and GTDB-Tk reference genomes. The classify workflow is recommended for obtaining taxonomic classifications, and this workflow only recommended if a *de novo* tree is desired. This workflow consists of five steps: *identify*, *align*, *infer*, *root*, and *decorate*. The *identify* and *align* steps are the same as in the classify workflow. The *infer* step uses [FastTree](http://www.microbesonline.org/fasttree/) with the WAG+GAMMA models to calculate a *de novo* tree. This tree is then rooted using a user specified outgroup and decorated with the GTDB taxonomy. 

The *de novo* workflow can be run as follows:
```
> gtdbtk de_novo_wf --genome_dir <my_genomes> --<marker_set> --outgroup_taxon <outgroup> --out_dir <output_dir>
```
This will process all genomes in <my_genomes> using the specified marker set and place the results in <output_dir>. Only genomes previously identified as being bacterial (archaeal) should be included when using the bacterial (archaeal) marker set. The tree will be rooted with the <outgroup> taxon. Identical to the classify workflow, the location of genomes can also be specified using a batch file with the --batchfile flag.

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

* Chaumeil PA, Hugenholtz P, Parks DH. 2018. GTDB-Tk: A toolkit to classify genomes with the Genome Taxonomy Database. \<in prep\>.

In the meantime, if you find the GTDB-Tk useful please cite this GitHub page. 

The GTDB taxonomy is described in:

* Parks DH, et al. 2018. [A standardized bacterial taxonomy based on genome phylogeny substantially revises the tree of life](https://www.nature.com/articles/nbt.4229). <i>Nat. Biotechnol.</i>, http://dx.doi.org/10.1038/nbt.4229
 
 We also strongly encourage you to cite the following 3rd party dependencies:

* Matsen FA, Kodner RB, Armbrust EV. 2010. [pplacer: linear time maximum-likelihood and Bayesian phylogenetic placement of sequences onto a fixed reference tree](https://www.ncbi.nlm.nih.gov/pubmed/21034504). <i>BMC Bioinformatics</i>, 11:538.
* Jain C, et al. 2017. [High-throughput ANI Analysis of 90K Prokaryotic Genomes Reveals Clear Species Boundaries](https://www.biorxiv.org/content/early/2017/11/27/225342). <i>bioRxiv</i>, https://doi.org/10.1101/225342.
* Hyatt D, et al. 2010. [Prodigal: prokaryotic gene recognition and translation initiation site identification](https://www.ncbi.nlm.nih.gov/pubmed/20211023). <i>BMC Bioinformatics</i>, 11:119. doi: 10.1186/1471-2105-11-119.
* Price MN, Dehal PS, Arkin AP. [FastTree 2 – Approximately Maximum-Likelihood Trees for Large Alignments](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2835736/). <i>PLoS One</i>, 5, e9490.
* Eddy SR. 2011. [Accelerated profile HMM searches](https://www.ncbi.nlm.nih.gov/pubmed/22039361). <i>PLOS Comp. Biol.</i>, 7:e1002195.

## Copyright

Copyright © 2017 Pierre-Alain Chaumeil. See LICENSE for further details.
