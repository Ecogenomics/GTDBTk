---
layout: default
title: Running GTDB-Tk
nav_order: 2
has_children: true
has_toc: false
---

# Running GTDB-Tk
{: .no_toc }

GTDB-Tk can be either installed locally, or run using the third-party web based platform KBase.



## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

## Hardware requirements
* ~100Gb of memory to run
* ~27Gb of storage
* ~1 hour per 1,000 genomes when using 64 CPUs

## Python libraries
GTDB-Tk is designed for Python >=3.6 and requires the following libraries which will be automatically installed:
* [dendropy](http://dendropy.org/)  >=4.1.0: a Python library for phylogenetic computing.
* [NumPy](https://numpy.org/) >=1.9.0: scientific computing with Python.


## Third-party software

GTDB-Tk makes use of the following 3rd party dependencies and assumes they are on your system path:
* [Prodigal](http://compbio.ornl.gov/prodigal/) >= 2.6.2: Hyatt D, et al. 2012. Gene and translation initiation site prediction in metagenomic sequences. <i>Bioinformatics</i>, 28, 2223-2230.
* [HMMER](http://hmmer.org/) >= 3.1: Eddy SR. 2011. Accelerated profile HMM searches. <i>PLoS Comp. Biol.</i>, 7, e1002195.
* [pplacer](http://matsen.fhcrc.org/pplacer/) >= 1.1: Matsen F, et al. 2010. pplacer: linear time maximum-likelihood and Bayesian phylogenetic placement of sequences onto a fixed reference tree. <i>BMC Bioinformatics</i>, 11, 538.
* [FastANI](https://github.com/ParBLiSS/FastANI) >= 1.2: Jain C, et al. 2018. High-throughput ANI analysis of 90K prokaryotic genomes reveals clear species boundaries. <i>Nature Communication</i>, 5114.
* [FastTree](http://www.microbesonline.org/fasttree/) >= 2.1.9: Price MN, et al. 2010 FastTree 2 -- Approximately Maximum-Likelihood Trees for Large Alignments. <i>PLoS ONE</i>, 5, e9490.
* [Mash](https://github.com/marbl/Mash) >= 2.2: Mash: fast genome and metagenome distance estimation using MinHash. Ondov BD, Treangen TJ, Melsted P, Mallonee AB, Bergman NH, Koren S, Phillippy AM. Genome Biol. 2016 Jun 20;17(1):132. doi: 10.1186/s13059-016-0997-x.

Please cite these tools if you use GTDB-Tk in your work.

## GTDB-Tk reference data

GTDB-Tk requires ~27G of external data that need to be downloaded and unarchived:
```
wget https://data.ace.uq.edu.au/public/gtdb/data/releases/release89/89.0/gtdbtk_r89_data.tar.gz
tar xvzf gtdbtk_r89_data.tar.gz
```

Reference data for prior releases of GTDB-Tk are available at:
```
wget https://data.ace.uq.edu.au/public/gtdbtk
```

