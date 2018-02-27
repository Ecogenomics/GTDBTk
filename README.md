# GTDB-Tk

[![version status](https://img.shields.io/pypi/v/gtdbtk.svg)](https://pypi.python.org/pypi/gtdbtk)
[![downloads](https://img.shields.io/pypi/dm/gtdbtk.svg)](https://pypi.python.org/pypi/gtdbtk)

GTDB-Tk is a software toolkit for assigning objective taxonomic classifications to bacterial and archaeal genomes. It is computationally 
efficient and designed to work with recent advances that allow hundreds or thousands of metagenome-assembled genomes (MAGs) to be obtained directly from environmental samples. However, it can equally be applied to isolate and single-cell genomes. The GTDB-Tk is open source and released under the GNU General Public License (Version 3).

Please note that GTDB-Tk is **still under active development** and is not ready for external users.
notifications about GTDB-Tk releases will be available through ACE Twitter account (https://twitter.com/ace_uq).

## Installation

GTDB-Tk requires 70G+ of disk.

GTDB-Tk requires the following Python libraries:
* [jinja2](http://jinja.pocoo.org/) >=2.7.3: a full featured template engine for Python.
* [mpld3](http://mpld3.github.io/) >= 0.2: D3 viewer for Matplotlib.
* [biolib](https://github.com/dparks1134/biolib) >= 0.0.44: Python package for common tasks in bioinformatic.
* [dendropy](http://dendropy.org/)  >= 4.1.0: A Python library for phylogenetics and phylogenetic computing: reading, writing, simulation, processing and manipulation of phylogenetic trees (phylogenies) and characters.
* [SciPy Stack](https://www.scipy.org/install.html): at least the Matplotlib, NumPy, and SciPy libraries

Jinja2, mpld3,dendropy and biolib should install as part of the GTDB-Tk if you install it via pip as described below. The SciPy Stack must be install seperately of GTDB-Tk.

GTDB-Tk makes use of the following 3rd party dependencies and assumes these are on your system path:
* [Prodigal](http://prodigal.ornl.gov/) >= 2.6.2: Hyatt D, et al. 2012. Gene and translation initiation site prediction in metagenomic sequences. <i>Bioinformatics</i>, 28, 2223-2230.
* [HMMER](http://http://hmmer.org/) >= 3.1: Eddy SR. 2011. Accelerated profile HMM searches. <i>PLoS Comp. Biol.</i>, 7, e1002195.
* [pplacer](http://matsen.fhcrc.org/pplacer/) >= 1.1: Matsen F, et al. 2010. pplacer: linear time maximum-likelihood and Bayesian phylogenetic placement of sequences onto a fixed reference tree. <i>BMC Bioinformatics</i>, 11, 538.
* [FastANI](https://github.com/ParBLiSS/FastANI) >= 1.0: Jain C, et al. 2017. High-throughput ANI Analysis of 90K Prokaryotic Genomes Reveals Clear Species Boundaries.<i>bioRxiv.</i> 256800.
* [FastTree](http://www.microbesonline.org/fasttree/) >= 2.1.9: Price MN, et al. 2010 FastTree 2 -- Approximately Maximum-Likelihood Trees for Large Alignments. <i>PLoS ONE</i>, 5, e9490.

GTDB-Tk also assumes the Python 2.7.x and Perl interpreters are on your system path.

GTDB-Tk requires external data that need to be downloaded and unarchived (preferably in the same folder):
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

GTDB-Tk requires a config file.
In the python lib/site-packages folder, go to the gtdbtk directory; once there :
```
cd config
cp config_template.py config.py
```
Edit the config.py file and modify different variables :
-GENERIC_PATH should point to the folder where data from the https://data.ace.uq.edu.au/public/gtdbtk/ has been downloaded. Make sure the variable finishes with a slash '/'.

## Quick Start

The functionality provided by GTDB-Tk can be accessed through the help menu:
```
> gtdbtk -h
```

Usage information about specific methods can also be accessed through the help menu, e.g.:
```
> gtdbtk classify_wf -h
```

## Classify Workflow

The classify workflow consists of three steps: identify, align, and classify. The identify step calls genes using [Prodigal](http://prodigal.ornl.gov/) and then uses HMM models and the [HMMER](http://http://hmmer.org/) package to identify the marker genes used for phylogenetic inference. As part of this search marker genes are aligned to their respective HMM model. The align step concatenates the aligned marker genes and applies all necessary filtering to the concatenated multiple sequence alignment. Finally, the classify step uses [pplacer](http://matsen.fhcrc.org/pplacer/) to find the maximum-likelihood placement of each genome into the GTDB-Tk reference tree based on its concatenated multiple sequence alignment. GTDB-Tk classify each genome based on its placement in the reference tree, relative evolutionary distance, and Mash distance (see Chaumeil PA et al., 2017 for details).
 
The classify workflow can be run as follows:
```
> gtdbtk classify_wf --genome_dir <my_genomes> --out_dir <output_dir>
```
This will process all genomes in <my_genomes> using both bacterial and archaeal marker sets and place the results in <output_dir>. Genomes must be in FASTA format. The location of genomes can also be specified using a batch file with the --batchfile flag. The batch file is simply a two column file indicating the location of each genome and the desired genome identified (i.e., a Newick compatible alphanumeric string). These fields must be seperated by a tab.

The workflow supports several optional flags, including:
* cpus: maximum number of CPUs to use
* proteins: indicates genome  files contain called proteins as amino acids and gene calling can be skipped

For other flags please consult the command line interface.

Here is an example run of this workflow:
```
> gtdbtk classify_wf --cpus 24 --genome_dir ./my_genomes --out_dir de_novo_output
```

##### Output files 
Each step of the Classify workflow generates a number of different files that would give users additional information on GTDB-Tk pipeline

Identify step:
* gtdbtk_bac120_markers_summary.tsv: Summary of unique markers, missing markers and unique markers of the 120 bacterial markers for each a submitted genomes.
* gtdbtk_ar122_markers_summary.tsv: Similar to gtdbtk_bac120_markers_summary.tsv but bases on the 122 archaeal markers.
* marker_genes folder: lists individual genome results for Gene calling unsing Prodigal and Gene matching based on TigrFam and Pfam markers.

Align step:
* *.user_msa.fasta: Multi sequence alignement  listing ONLY the submitted genomes
* *.msa.fasta: Multi sequence alignement  listing ALL genomes ( submitted + reference genomes)
* *.filtered.tsv: Display the list of genomes with an insufficient number of amino acids in MSA.

Classify step:
* *.classification_pplacer.tsv: Classification of user genomes based only on pplacer.
* *.classification.tsv: Classification of user genomes based on the FastANI, RED Value and pplacer. 
* *.fastani_results.tsv: Tab deleimited file listing the user genome, the reference genome and the Average nucleotide identity between the two
* *.classify.tree: Reference tree with all user genomes placed with pplacer
* *.summary.tsv: if user genomes have been classified using FastANI,Red values or only the topology of the tree.
* *.red_dictionary: Median Red values for the ranks Phylum to Genus




## Validating Species Assignments

The GTDB-Tk uses FastANI to estimate the ANI between genomes.We recommend that species assignments made by the GTDB-Tk be validate using ANI distances against a suitably set of reference genomes. ANI values can be calculated using [ANIcalculator](https://ani.jgi-psf.org/html/home.php).

## De Novo Workflow

The <i>de novo</i> workflow infers a new tree containing all reference GTDB-Tk genomes and user supply genomes. The classify workflow is recommended for obtaining taxonomic classifications, and this workflow only recommended if a <i>de novo</i> tree is desired. This workflow consists of five steps: identify, align, infer, root, and decorate. The first to are the same as in the classify workflow. The infer step uses [FastTree](http://www.microbesonline.org/fasttree/) with the WAG+GAMMA models to calculate a <i>de novo</i> tree. This tree is then rooted using a user specified outgroup and then decorated with the GTDB taxonomy. 

The classify workflow can be run as follows:
```
> gtdbtk de_novo_wf --genome_dir <my_genomes> --<marker_set> --outgroup_taxon <outgroup> --out_dir <output_dir>
```
This will process all genomes in <my_genomes> using the specified marker set and place the results in <output_dir>. The tree will thbe rooted with the <outgroup> taxon. Identical to the classify workflow, the location of genomes can also be specified using a batch file with the --batchfile flag and any of the 3 marker sets defined by GTDB-Tk used to infer the tree.

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

If you find this package useful, please cite:

Chaumeil PA, Parks DH, Hugenholtz P. 2017. GTDB-Tk: A toolkit to classify genomes with the Genome Taxonomy Database. <in prep>.


## Copyright

Copyright © 2017 Pierre-Alain Chaumeil. See LICENSE for further details.
