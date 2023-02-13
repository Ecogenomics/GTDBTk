# GTDB-Tk

[![PyPI](https://img.shields.io/pypi/v/gtdbtk.svg)](https://pypi.python.org/pypi/gtdbtk)
[![PyPI Downloads](https://pepy.tech/badge/gtdbtk)](https://pepy.tech/project/gtdbtk)
[![Bioconda](https://img.shields.io/conda/vn/bioconda/gtdbtk.svg?color=43b02a)](https://anaconda.org/bioconda/gtdbtk)
[![BioConda Downloads](https://img.shields.io/conda/dn/bioconda/gtdbtk.svg?style=flag&label=downloads&color=43b02a)](https://anaconda.org/bioconda/gtdbtk)
[![Docker Image Version (latest by date)](https://img.shields.io/docker/v/ecogenomic/gtdbtk?sort=date&color=299bec&label=docker)](https://hub.docker.com/r/ecogenomic/gtdbtk)
[![Docker Pulls](https://img.shields.io/docker/pulls/ecogenomic/gtdbtk?color=299bec&label=pulls)](https://hub.docker.com/r/ecogenomic/gtdbtk)

GTDB-Tk is a software toolkit for assigning objective taxonomic classifications to bacterial and archaeal genomes based 
on the Genome Database Taxonomy ([GTDB](https://gtdb.ecogenomic.org/)). It is designed to work with recent advances that 
allow hundreds or thousands of metagenome-assembled genomes (MAGs) to be obtained directly from environmental samples. 
It can also be applied to isolate and single-cell genomes. The GTDB-Tk is open source and released under the 
[GNU General Public License (Version 3)](https://www.gnu.org/licenses/gpl-3.0.en.html).

Notifications about GTDB-Tk releases will be available through the [GTDB Twitter](https://twitter.com/ace_gtdb) 
account and the [GTDB Announcements Forum](https://forum.gtdb.ecogenomic.org/c/announcements/10).

Please post questions and issues related to GTDB-Tk on the Issues section of the GitHub repository. Questions 
related to the [GTDB](https://gtdb.ecogenomic.org/) can be posted on the [GTDB Forum](https://forum.gtdb.ecogenomic.org/) 
or sent to the [GTDB team](https://gtdb.ecogenomic.org/about).


## 🚀 Getting started

Be sure to check the [hardware requirements](https://ecogenomics.github.io/GTDBTk/installing/index.html), then choose your preferred method:

* [Bioconda](https://ecogenomics.github.io/GTDBTk/installing/bioconda.html)
* [Docker](https://ecogenomics.github.io/GTDBTk/installing/docker.html)
* [pip](https://ecogenomics.github.io/GTDBTk/installing/pip.html)


## 📖 Documentation

Documentation for GTDB-Tk can be found [here](https://ecogenomics.github.io/GTDBTk/).


## ✨ New Features

GTDB-Tk v2.2.0+ includes the following new features:
- GTDB-TK `classify` and `classify_wf` have changed in version 2.2.0+. There is now an ANI classification stage (`ANI screen`) that precedes classification by placement in a reference tree.
  - **This is now the default behavior for `classify` and `classify_wf`.**
  - In `classify`, user genomes are first compared against a Mash database comprised of all GTDB representative genomes and genome pairs of sufficient similarity processed by FastANI. User genomes classified to a GTDB representative based on FastANI results are not run through pplacer. 
  - In the `classify_wf` workflow, genomes are classified using Mash and FastANI before executing the identify step. User genomes classified with FastANI are not run through the remainder of the pipeline (identify, align, classify).
  - To classify genomes without the additional `ani_screen` step, use the `--skip_ani_screen` flag.

## 📈 Performance
Using ANI screen "can" reduce computation by >50%, although it depends on the set of input genomes. A set of input genomes consisting primarily of new species will not benefit from ANI screen as much as a set of genomes that are largely assigned to GTDB species clusters. In the latter case, the ANI screen will reduce the number of genomes that need to be classified by pplacer which reduces computation time subsantially (between 25% and 60% in our testing).

## 📚 References

GTDB-Tk is described in:

* Chaumeil PA, et al. 2022. [GTDB-Tk v2: memory friendly classification with the Genome Taxonomy Database](https://academic.oup.com/bioinformatics/advance-article-abstract/doi/10.1093/bioinformatics/btac672/6758240?utm_source=advanceaccess&utm_campaign=bioinformatics&utm_medium=email). <i>Bioinformatics</i>, btac672.
* Chaumeil PA, et al. 2019. [GTDB-Tk: A toolkit to classify genomes with the Genome Taxonomy Database](https://academic.oup.com/bioinformatics/advance-article-abstract/doi/10.1093/bioinformatics/btz848/5626182). <i>Bioinformatics</i>, btz848.

The Genome Taxonomy Database (GTDB) is described in:

* Parks, D.H., et al. (2021). [GTDB: an ongoing census of bacterial and archaeal diversity through a phylogenetically consistent, rank normalized and complete genome-based taxonomy](https://academic.oup.com/nar/advance-article/doi/10.1093/nar/gkab776/6370255). <i>Nucleic Acids Research</i>, <b>50</b>: D785–D794.
* Rinke, C, et al. (2021). [A standardized archaeal taxonomy for the Genome Taxonomy Database](https://www.nature.com/articles/s41564-021-00918-8). <i>Nature Microbiology</i>, <b>6</b>: 946–959.
* Parks, D.H., et al. 2020. [A complete domain-to-species taxonomy for Bacteria and Archaea](https://rdcu.be/b3OI7). <i>Nature Biotechnology</i>, https://doi.org/10.1038/s41587-020-0501-8.
* Parks DH, et al. 2018. [A standardized bacterial taxonomy based on genome phylogeny substantially revises the tree of life](https://www.nature.com/articles/nbt.4229). <i>Nature Biotechnology</i>, http://dx.doi.org/10.1038/nbt.4229.
 

We strongly encourage you to cite the following 3rd party dependencies:

* Matsen FA, et al. 2010. [pplacer: linear time maximum-likelihood and Bayesian phylogenetic placement of sequences onto a fixed reference tree](https://www.ncbi.nlm.nih.gov/pubmed/21034504). <i>BMC Bioinformatics</i>, 11:538.
* Jain C, et al. 2019. [High-throughput ANI Analysis of 90K Prokaryotic Genomes Reveals Clear Species Boundaries](https://www.nature.com/articles/s41467-018-07641-9). <i>Nat. Communications</i>, doi: 10.1038/s41467-018-07641-9.
* Hyatt D, et al. 2010. [Prodigal: prokaryotic gene recognition and translation initiation site identification](https://www.ncbi.nlm.nih.gov/pubmed/20211023). <i>BMC Bioinformatics</i>, 11:119. doi: 10.1186/1471-2105-11-119.
* Price MN, et al. 2010. [FastTree 2 - Approximately Maximum-Likelihood Trees for Large Alignments](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2835736/). <i>PLoS One</i>, 5, e9490.
* Eddy SR. 2011. [Accelerated profile HMM searches](https://www.ncbi.nlm.nih.gov/pubmed/22039361). <i>PLOS Comp. Biol.</i>, 7:e1002195.
* Ondov BD, et al. 2016. [Mash: fast genome and metagenome distance estimation using MinHash](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0997-x). <i>Genome Biol</i> 17, 132. doi: 10.1186/s13059-016-0997-x.


## © Copyright

Copyright 2017 Pierre-Alain Chaumeil. See LICENSE for further details.
