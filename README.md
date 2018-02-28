# GTDB-Tk

[![version status](https://img.shields.io/pypi/v/gtdbtk.svg)](https://pypi.python.org/pypi/gtdbtk)

GTDB-Tk is a software toolkit for assigning objective taxonomic classifications to bacterial and archaeal genomes. It is computationally 
efficient and designed to work with recent advances that allow hundreds or thousands of metagenome-assembled genomes (MAGs) to be obtained directly from environmental samples. It can also be applied to isolate and single-cell genomes. The GTDB-Tk is open source and released under the GNU General Public License (Version 3).

GTDB-Tk is **under active development and validation**. Please independently confirm the GTDB-Tk predictions by manually inspecting the tree and bringing any discrepencies to our attention. Notifications about GTDB-Tk releases will be available through the ACE Twitter account (https://twitter.com/ace_uq).

## Installation

GTDB-Tk requires the following Python libraries:
* [jinja2](http://jinja.pocoo.org/) >=2.7.3: a full featured template engine for Python.
* [mpld3](http://mpld3.github.io/) >= 0.2: D3 viewer for Matplotlib.
* [biolib](https://github.com/dparks1134/biolib) >= 0.0.44: Python package for common tasks in bioinformatic.
* [dendropy](http://dendropy.org/)  >= 4.1.0: A Python library for phylogenetics and phylogenetic computing: reading, writing, simulation, processing and manipulation of phylogenetic trees (phylogenies) and characters.
* [SciPy Stack](https://www.scipy.org/install.html): at least the Matplotlib, NumPy, and SciPy libraries

Jinja2, mpld3, dendropy and biolib will be install as part of GTDB-Tk when installing via pip as described below. The SciPy Stack must be install seperately.

GTDB-Tk makes use of the following 3rd party dependencies and assumes these are on your system path:
* [Prodigal](http://prodigal.ornl.gov/) >= 2.6.2: Hyatt D, et al. 2012. Gene and translation initiation site prediction in metagenomic sequences. <i>Bioinformatics</i>, 28, 2223-2230.
* [HMMER](http://http://hmmer.org/) >= 3.1: Eddy SR. 2011. Accelerated profile HMM searches. <i>PLoS Comp. Biol.</i>, 7, e1002195.
* [pplacer](http://matsen.fhcrc.org/pplacer/) >= 1.1: Matsen F, et al. 2010. pplacer: linear time maximum-likelihood and Bayesian phylogenetic placement of sequences onto a fixed reference tree. <i>BMC Bioinformatics</i>, 11, 538.
* [FastANI](https://github.com/ParBLiSS/FastANI) >= 1.0: Jain C, et al. 2017. High-throughput ANI Analysis of 90K Prokaryotic Genomes Reveals Clear Species Boundaries.<i>bioRxiv.</i> 256800.
* [FastTree](http://www.microbesonline.org/fasttree/) >= 2.1.9: Price MN, et al. 2010 FastTree 2 -- Approximately Maximum-Likelihood Trees for Large Alignments. <i>PLoS ONE</i>, 5, e9490.

GTDB-Tk also assumes the Python 2.7.x and Perl interpreters are on your system path.

GTDB-Tk requires ~70G+ of external data that need to be downloaded and unarchived (preferably in the same directory):
```
wget https://data.ace.uq.edu.au/public/gtdbtk/release_80/fastani.tar.gz
wget https://data.ace.uq.edu.au/public/gtdbtk/release_80/markers.tar.gz
wget https://data.ace.uq.edu.au/public/gtdbtk/release_80/masks.tar.gz
wget https://data.ace.uq.edu.au/public/gtdbtk/release_80/msa.tar.gz
wget https://data.ace.uq.edu.au/public/gtdbtk/release_80/pplacer.tar.gz
wget https://data.ace.uq.edu.au/public/gtdbtk/release_80/taxonomy.tar.gz
```

Once these are installed, GTDB-Tk can be installed using [pip](https://pypi.python.org/pypi/gtdbtk):
```
> pip install gtdbtk
```

GTDB-Tk requires a config file. In the Python lib/site-packages directory, go to the gtdbtk directory and setup this config file:
```
cd config
cp config_template.py config.py
```
Edit the config.py file and modify different variables:
-GENERIC_PATH should point to the directory containing the data downloaded from the https://data.ace.uq.edu.au/public/gtdbtk/. Make sure the variable finishes with a slash '/'.

## Quick Start

The functionality provided by GTDB-Tk can be accessed through the help menu:
```
> gtdbtk -h
```

Usage information about each methods can also be accessed through their species help menu, e.g.:
```
> gtdbtk classify_wf -h
```

## Classify Workflow

The classify workflow consists of three steps: *identify*, *align*, and *classify*. The *identify* step calls genes using [Prodigal](http://prodigal.ornl.gov/) and then uses HMM models and the [HMMER](http://http://hmmer.org/) package to identify the marker genes used for phylogenetic inference. Consistent alignments are obtained by aligning marker genes to their respective HMM model. The *align* step concatenates the aligned marker genes and applies all necessary filtering to the concatenated multiple sequence alignment (MSA). Finally, the *classify* step uses [pplacer](http://matsen.fhcrc.org/pplacer/) to find the maximum-likelihood placement of each genome's concatenated protein alignment in the GTDB-Tk reference tree. GTDB-Tk classifies each genome based on its placement in the reference tree, its relative evolutionary distance, and FastANI distance (see Chaumeil PA et al., 2018 for details).
 
The classify workflow can be run as follows:
```
> gtdbtk classify_wf --genome_dir <my_genomes> --out_dir <output_dir>
```
This will process all genomes in <my_genomes> using both bacterial and archaeal marker sets and place the results in <output_dir>. Genomes must be in FASTA format. The location of genomes can also be specified using a batch file with the --batchfile flag. The batch file is simply a two column file indicating the location of each genome and the desired genome identifier (i.e., a Newick compatible alphanumeric string). These fields must be seperated by a tab.

The workflow supports several optional flags, including:
* cpus: maximum number of CPUs to use

For other flags please consult the command line interface.

Here is an example run of this workflow:
```
> gtdbtk classify_wf --cpus 24 --genome_dir ./my_genomes --out_dir gtdbtk_output
```

The taxonomic classification of each bacterial and archaeal genome is contained in the \<prefix\>.bac120.classification.tsv and \<prefix\>.ar122.classification.tsv output files.

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
* \<prefix\>.classification.tsv: classification of user genomes based on the FastANI, RED values, and pplacer. This is the primary output of the GTDB-Tk and contains the taxonomic classification we recommend
 * \<prefix\>.summary.tsv: additional information regarding the criteria used to classify a genome
* \<prefix\>.classification_pplacer.tsv: classification of user genomes based only on pplacer
* \<prefix\>.fastani_results.tsv: tab delimited file indicating the average nucleotide identity (ANI) between a user and reference genome for cases where classification is based on the ANI statistic
* \<prefix\>.classify.tree: reference tree in Newick format containing all user genomes placed with pplacer in the GTDB-Tk reference tree
* \<prefix\>.red_dictionary: median RED values for taxonomic ranks

## Validating Species Assignments

The GTDB-Tk uses FastANI to estimate the average nucleotide identity (ANI) between genomes. Species assignments are made using an ANI criteria of 95%. Information about species assignments can be found in the <prefix>.fastani_results.tsv output file.

## De Novo Workflow
**under active development**
The *de novo* workflow infers a new tree containing all user supply and GTDB-Tk reference genomes. The classify workflow is recommended for obtaining taxonomic classifications, and this workflow only recommended if a *de novo* tree is desired. This workflow consists of five steps: *identify*, *align*, *infer*, *root*, and *decorate*. The *identify* and *align* steps are the same as in the classify workflow. The *infer* step uses [FastTree](http://www.microbesonline.org/fasttree/) with the WAG+GAMMA models to calculate a *de novo* tree. This tree is then rooted using a user specified outgroup and decorated with the GTDB taxonomy. 

The *de novo* workflow can be run as follows:
```
> gtdbtk de_novo_wf --genome_dir <my_genomes> --<marker_set> --outgroup_taxon <outgroup> --out_dir <output_dir>
```
This will process all genomes in <my_genomes> using the specified marker set and place the results in <output_dir>. Only genomes previously identified as being bacterial (archaeal) should be included when using the bacterial (archaeal) marker set. The tree will be rooted with the <outgroup> taxon. Identical to the classify workflow, the location of genomes can also be specified using a batch file with the --batchfile flag.

The workflow supports several optional flags, including:
* cpus: maximum number of CPUs to use
* min_perc_aa: filter genomes with an insufficient percentage of AA in the MSA (default: 50)
* taxa_filter: filter genomes to taxa within specific taxonomic groups
* prot_model:  protein substitution model for tree inference (LG or WAG; default: WAG)

For other flags please consult the command line interface.

Here is an example run of this workflow:
```
> gtdbtk de_novo_wf --genome_dir ./genomes --bac120_ms --outgroup_taxon p__Acetothermia --taxa_filter p__Firmicutes --out_dir de_novo_output
```

## Individual Steps

All steps comprising the classify and <i>de novo</i> workflows can be run independently if desired. Please consult the command line interface for specific details on running each of these steps.

## Cite

A manuscript describing the GTDB-Tk is currently being prepared:

Chaumeil PA, Hugenholtz P, Parks DH. 2018. GTDB-Tk: A toolkit to classify genomes with the Genome Taxonomy Database. \<in prep\>.
 
In the meantime, if you find the GTDB-Tk useful please cite this GitHub page. Please also consider citing the 3rd party applications required by GTDB-Tk such as [Prodigal](http://prodigal.ornl.gov/), [HMMER](http://http://hmmer.org/), [pplacer](http://matsen.fhcrc.org/pplacer/), [FastANI](https://github.com/ParBLiSS/FastANI), and
[FastTree](http://www.microbesonline.org/fasttree/).

## Copyright

Copyright © 2017 Pierre-Alain Chaumeil. See LICENSE for further details.
