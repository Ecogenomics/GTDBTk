# GTDB-Tk

[![PyPI](https://img.shields.io/pypi/v/gtdbtk.svg)](https://pypi.python.org/pypi/gtdbtk)
[![PyPI Downloads](https://pepy.tech/badge/gtdbtk)](https://pepy.tech/project/gtdbtk)
[![Bioconda](https://img.shields.io/conda/vn/bioconda/gtdbtk.svg?color=43b02a)](https://anaconda.org/bioconda/gtdbtk)
[![BioConda Downloads](https://img.shields.io/conda/dn/bioconda/gtdbtk.svg?style=flag&label=downloads&color=43b02a)](https://anaconda.org/bioconda/gtdbtk)
[![Docker Image Version (latest by date)](https://img.shields.io/docker/v/ecogenomic/gtdbtk?sort=date&color=299bec&label=docker)](https://hub.docker.com/r/ecogenomic/gtdbtk)
[![Docker Pulls](https://img.shields.io/docker/pulls/ecogenomic/gtdbtk?color=299bec&label=pulls)](https://hub.docker.com/r/ecogenomic/gtdbtk)

<b>[GTDB-Tk v2.1.0](https://ecogenomics.github.io/GTDBTk/announcements.html) was released on May 11, 2022. Upgrading is recommended.</b>  
<b> Please note v2.1.0+ is not compatible with GTDB-Tk package [R207_v1](https://data.gtdb.ecogenomic.org/releases/latest/auxillary_files/gtdbtk_data.tar.gz). It is necessary to upgrade to GTDB-Tk package [R207_v2](https://data.gtdb.ecogenomic.org/releases/latest/auxillary_files/gtdbtk_v2_data.tar.gz).</b>

GTDB-Tk is a software toolkit for assigning objective taxonomic classifications to bacterial and archaeal genomes based on the Genome Database Taxonomy ([GTDB](https://gtdb.ecogenomic.org/)). It is designed to work with recent advances that allow hundreds or thousands of metagenome-assembled genomes (MAGs) to be obtained directly from environmental samples. It can also be applied to isolate and single-cell genomes. The GTDB-Tk is open source and released under the [GNU General Public License (Version 3)](https://www.gnu.org/licenses/gpl-3.0.en.html).

Notifications about GTDB-Tk releases will be available through the [GTDB Twitter](https://twitter.com/ace_gtdb) account and the [GTDB Announcements Forum](https://forum.gtdb.ecogenomic.org/c/announcements/10).

Please post questions and issues related to GTDB-Tk on the Issues section of the GitHub repository. Questions related to the [GTDB](https://gtdb.ecogenomic.org/) can be posted on the [GTDB Forum](https://forum.gtdb.ecogenomic.org/) or sent to the [GTDB team](https://gtdb.ecogenomic.org/about).

## New Features

GTDB-Tk v2.1.0 includes the following new features:
- GTDB-TK now uses a **divide-and-conquer** approach where the bacterial reference tree is split into multiple **class**-level subtrees. This reduces the memory requirements of GTDB-Tk from **320 GB** of RAM when using the full GTDB R07-RS207 reference tree to approximately **55 GB**. A manuscript describing this approach is in preparation. If you wish to continue using the full GTDB reference tree use the `--full-tree` flag.  
This is the main change from v2.0.0. The split tree approach has been modified from order-level trees to class-level trees to resolve specific classification issues (See [#383](https://github.com/Ecogenomics/GTDBTk/issues/383)). 
- Genomes that cannot be assigned to a domain (e.g. genomes with no bacterial or archaeal markers or genomes with no genes called by Prodigal) are now reported in the `gtdbtk.bac120.summary.tsv` as 'Unclassified'
- Genomes filtered out during the alignment step are now reported in the `gtdbtk.bac120.summary.tsv` or `gtdbtk.ar53.summary.tsv` as 'Unclassified Bacteria/Archaea'
- `--write_single_copy_genes` flag in now available in the `classify_wf` and `de_novo_wf` workflows.


## Documentation
Documentation for GTDB-Tk can be found [here](https://ecogenomics.github.io/GTDBTk/).

## References

GTDB-Tk is described in:

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

## Copyright

Copyright 2017 Pierre-Alain Chaumeil. See LICENSE for further details.
