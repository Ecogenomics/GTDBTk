# GTDB-Tk

[![PyPI](https://img.shields.io/pypi/v/gtdbtk.svg)](https://pypi.python.org/pypi/gtdbtk)
[![PyPI Downloads](https://pepy.tech/badge/gtdbtk)](https://pepy.tech/project/gtdbtk)
[![Bioconda](https://img.shields.io/conda/vn/bioconda/gtdbtk.svg?color=43b02a)](https://anaconda.org/bioconda/gtdbtk)
[![BioConda Downloads](https://img.shields.io/conda/dn/bioconda/gtdbtk.svg?style=flag&label=downloads&color=43b02a)](https://anaconda.org/bioconda/gtdbtk)
[![Docker Image Version (latest by date)](https://img.shields.io/docker/v/ecogenomic/gtdbtk?sort=date&color=299bec&label=docker)](https://hub.docker.com/r/ecogenomic/gtdbtk)
[![Docker Pulls](https://img.shields.io/docker/pulls/ecogenomic/gtdbtk?color=299bec&label=pulls)](https://hub.docker.com/r/ecogenomic/gtdbtk)

<b>[GTDB-Tk v2.0.0](https://ecogenomics.github.io/GTDBTk/announcements.html) was released on April 8, 2022 along with new reference data for [GTDB R07-RS207](https://gtdb.ecogenomic.org/). Upgrading is recommended.</b>  
<b> Please note v2.0.0+ is not compatible with GTDB R06-RS202. </b>

GTDB-Tk is a software toolkit for assigning objective taxonomic classifications to bacterial and archaeal genomes based on the Genome Database Taxonomy [GTDB](https://gtdb.ecogenomic.org/). It is designed to work with recent advances that allow hundreds or thousands of metagenome-assembled genomes (MAGs) to be obtained directly from environmental samples. It can also be applied to isolate and single-cell genomes. The GTDB-Tk is open source and released under the [GNU General Public License (Version 3)](https://www.gnu.org/licenses/gpl-3.0.en.html).

Notifications about GTDB-Tk releases will be available through the GTDB Twitter account (https://twitter.com/ace_gtdb).

Please post questions and issues related to GTDB-Tk on the Issues section of the GitHub repository. Questions related to the [GTDB](https://gtdb.ecogenomic.org/) can be posted on the [GTDB Forum](https://forum.gtdb.ecogenomic.org/) or sent to the [GTDB team](https://gtdb.ecogenomic.org/about).

## New Features

GTDB-Tk v2.0.0 includes the following new features:
- GTDB-TK now uses a **divide-and-conquer** approach where the bacterial reference tree is split into multiple order-level subtrees. This reduces the memory requirements of GTDB-Tk from **320 GB** of RAM when using the full GTDB R07-RS207 reference tree to approximately **35 GB**. A manuscript describing this approach is in preparation. If you wish to continue using the full GTDB reference tree use the `--full-tree` flag.
- Archaeal classification now uses a refined set of 53 archaeal-specific marker genes based on the recent publication by [Dombrowski et al., 2020](https://www.nature.com/articles/s41467-020-17408-w). This set of archaeal marker genes is now used by GTDB for curating the archaeal taxonomy.
- By default, all directories containing intermediate results are **now removed** by default at the end of the `classify_wf` and `de_novo_wf` pipelines. If you wish to retain these intermediates files use the `--keep-intermediates` flag.
- All MSA files produced by the `align` step are now compressed with gzip.
- The classification summary and failed genomes files are now the only files linked in the root directory of `classify_wf`.


## Documentation
Documentation for GTDB-Tk can be found [here](https://ecogenomics.github.io/GTDBTk/).

## References

GTDB-Tk is described in:

* Chaumeil PA, et al. 2019. [GTDB-Tk: A toolkit to classify genomes with the Genome Taxonomy Database](https://academic.oup.com/bioinformatics/advance-article-abstract/doi/10.1093/bioinformatics/btz848/5626182). <i>Bioinformatics</i>, btz848.

The Genome Taxonomy Database (GTDB) is described in:

* Parks, D.H., et al. 2020. [A complete domain-to-species taxonomy for Bacteria and Archaea](https://rdcu.be/b3OI7). <i>Nature Biotechnology</i>, https://doi.org/10.1038/s41587-020-0501-8.

* Parks DH, et al. 2018. [A standardized bacterial taxonomy based on genome phylogeny substantially revises the tree of life](https://www.nature.com/articles/nbt.4229). <i>Nature Biotechnology</i>, http://dx.doi.org/10.1038/nbt.4229.
 

We strongly encourage you to cite the following 3rd party dependencies:

* Matsen FA, et al. 2010. [pplacer: linear time maximum-likelihood and Bayesian phylogenetic placement of sequences onto a fixed reference tree](https://www.ncbi.nlm.nih.gov/pubmed/21034504). <i>BMC Bioinformatics</i>, 11:538.
* Jain C, et al. 2019. [High-throughput ANI Analysis of 90K Prokaryotic Genomes Reveals Clear Species Boundaries](https://www.nature.com/articles/s41467-018-07641-9). <i>Nat. Communications</i>, doi: 10.1038/s41467-018-07641-9.
* Hyatt D, et al. 2010. [Prodigal: prokaryotic gene recognition and translation initiation site identification](https://www.ncbi.nlm.nih.gov/pubmed/20211023). <i>BMC Bioinformatics</i>, 11:119. doi: 10.1186/1471-2105-11-119.
* Price MN, et al. 2010. [FastTree 2 - Approximately Maximum-Likelihood Trees for Large Alignments](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2835736/). <i>PLoS One</i>, 5, e9490.
* Eddy SR. 2011. [Accelerated profile HMM searches](https://www.ncbi.nlm.nih.gov/pubmed/22039361). <i>PLOS Comp. Biol.</i>, 7:e1002195.
* Ondov BD, et al. 2016. [Mash: fast genome and metagenome distance estimation using MinHash](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0997-x). <i>Genome Biol</i> 17, 132. doi: doi: 10.1186/s13059-016-0997-x.

## Copyright

Copyright 2017 Pierre-Alain Chaumeil. See LICENSE for further details.
