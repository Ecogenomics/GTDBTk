# GTDB-Tk

[![version status](https://img.shields.io/pypi/v/gtdbtk.svg)](https://pypi.python.org/pypi/gtdbtk)
[![downloads](https://img.shields.io/pypi/dm/gtdbtk.svg)](https://pypi.python.org/pypi/gtdbtk)

GTDB-Tk is a software toolkit for assigning objective taxonomic classifications to bacterial and archaeal genomes. It is computationally 
efficient and designed to work with recent advances that allow hundreds or thousands of metagenome-assembled genomes (MAGs) to be obtained 
directly from environmental samples. However, it can equally be applied to isolate and single-cell genomes. The GTDB-Tk is open source and 
released under the GNU General Public License (Version 3).

## Installation

GTDB-Tk makes use of the following 3rd party dependencies and assumes these are on your system path:
* [prodigal](http://prodigal.ornl.gov/) >= 2.6.2: Hyatt D, et al. 2012. Gene and translation initiation site prediction in metagenomic sequences. <i>Bioinformatics</i> 28: 2223-2230.
* [hmmer](http://http://hmmer.org/) >= 3.1: Eddy SR. 2011. Accelerated profile HMM searches. PLoS Comp. Biol, 7, e1002195.
* [pplacer](http://matsen.fhcrc.org/pplacer/) >= 1.1: Matsen, F. et al. 2010. pplacer: linear time maximum-likelihood and Bayesian phylogenetic placement of sequences onto a fixed reference tree. BMC Bioinformatics, 11, 538.

Once these are installed, GTDB-Tk can be installed using [pip](https://pypi.python.org/pypi/gtdbtk):
```
> sudo pip install gtdbtk
```

## Quick Start

The functionality provided by GTDB-Tk can be accessed through the help menu:
```
> gtdbtk -h
```

Usage information about specific methods can also be accessed through the help menu, e.g.:
```
> gtdbtk classify_wf –h
```

## Classify Workflow

## De novo workflow

## Validating Species Assignments

The GTDB-Tk uses Mash distances to estimate the ANI between genomes. Mash is computationally favourable, but is primarily used for practical considerations as it removes the need to have all references genomes avaliable to the GTDB-Tk in order to calculate ANI values. We recommend that species assignments made by the GTDB-Tk be validate using ANI distances against a suitably large set of reference genomes. ANI values can be calculated using [gANI](https://ani.jgi-psf.org/html/home.php).

## Cite

If you find this package useful, please cite:

<manuscript under preperation>


## Copyright

Copyright © 2017 Pierre-Alain Chaumeil. See LICENSE for further details.
