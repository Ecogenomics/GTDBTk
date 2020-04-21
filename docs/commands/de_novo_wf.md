---
layout: default
title: de_novo_wf
nav_order: 3
parent: Commands
---

# *de*_*novo*_wf 

In development
{: .label .label-yellow }

For arguments and output files, see each of the individual steps:

* [infer](infer.html)
* [align](align.html)
* [infer](infer.html)
* [root](root.html)
* [decorate](decorate.html)


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


## Example

### Input

```bash
gtdbtk de_novo_wf --genome_dir genomes/ --outgroup_taxon p__Nanobacteria --ar122_ms --out_dir de_novo_wf --cpus 3

gtdbtk de_novo_wf --genome_dir ./genomes --bac120_ms --outgroup_taxon p__Chloroflexota --taxa_filter p__Firmicutes --out_dir de_novo_output
```
