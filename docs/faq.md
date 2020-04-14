---
layout: default
title: FAQ
nav_order: 3
---


# FAQ

### Why is there a discrepancy in the naming system between GTDB-Tk and NCBI or Silva taxonomic names?  

GTDB-Tk uses the GTDB taxonomy ([http://gtdb.ecogenomic.org/](http://gtdb.ecogenomic.org/)). 
This taxonomy is similar, but not identical to NCBI and Silva. 
In many cases the GTDB taxonomy more strictly follows the nomenclatural rules for rank suffixes which is why there is Nitrospirota instead of Nitrospirae.

### GTDB-Tk crashes after reaching the job memory limit on my Job Scheduler / why does pplacer fail?    

We believe that this is an issue with how Linux reports memory usage when `pplacer` is run with multiple CPUs. 
It appears that Linux believes that each thread is consuming the total amount of memory allocated, not sharing it. We have typically oberseved this issue affecting queuing systems/HPCs.

A solution is to run the classify step using the `--pplacer_cpus` argument to reduce the number of threads allocated to pplacer.


## Validating species assignments with average nucleotide identity

The GTDB-Tk uses [FastANI](https://github.com/ParBLiSS/FastANI) to estimate the ANI between genomes. We  recommend you have FastANI >= 1.2 as this version introduces a fix that makes the results deterministic. A query genome is only classified as belonging to the same species as a reference genome if the ANI between the genomes is within the species ANI circumscription radius (typically, 95%) and the alignment fraction (AF) is >=0.65. In some circumstances, the phylogenetic placement of a query genome may not support the species assignment. GTDB r89 strictly uses ANI to circumscribe species and GTDB-Tk follows this methodology. The species-specific ANI circumscription radii are avaliable from the [GTDB](https://gtdb.ecogenomic.org/) website.