---
layout: default
title: align
nav_order: 5
parent: Commands
---

# align

## Arguments

### Required

| Argument     | Description                |
|:------------|:---------------------------|
| `--identify_dir IDENTIFY_DIR`          | output directory of `identify` command |
| `--out_dir OUT_DIR`          | directory to output files |


### Optional

| Argument     | Description                |
|:------------|:---------------------------|
| `--skip_gtdb_refs`          | do not include GTDB reference genomes in multiple sequence alignment (default: `False`) |
| `--taxa_filter TAXA_FILTER`          | filter GTDB genomes to taxa (comma separated) within specific taxonomic groups (e.g., `d__Bacteria` or `p__Proteobacteria, p__Actinobacteria`) |
| `--min_perc_aa MIN_PERC_AA`          | filter genomes with an insufficient percentage of AA in the MSA (inclusive bound) (default: `10`) |
| `--cols_per_gene COLS_PER_GENER`          | maximum number of columns to retain per gene (default: `42`) |
| `--min_consensus MIN_CONSENSUS`          | minimum percentage of the same amino acid required to retain column (inclusive bound) (default: `25`) |
| `--max_consensus MAX_CONSENSUS`          | maximum percentage of the same amino acid required to retain column (exclusive bound) (default: `95`) |
| `--min_perc_taxa MIN_PERC_TAXA`          | minimum percentage of taxa required to retain column (inclusive bound) (default: `50`) |
| `--rnd_seed RND_SEED`          | random seed to use for selecting columns |
| `--prefix PREFIX`          | desired prefix for output files (default: `gtdbtk`) |
| `--cpus CPUS`          | number of CPUs to use (default: `1`) |
| `--debug`          | create intermediate files for debugging purposes (default: `False`) |
| `-h`, `--help`          | show help message |


#### Mutually exclusive

| Argument     | Description                |
|:------------|:---------------------------|
| `--custom_msa_filters`          | perform custom filtering of MSA with cols_per_gene, min_consensus max_consensus, and min_perc_taxa parameters instead of using canonical mask (default: `False`) |
| `--skip_trimming`          | skip trimming step and return the full MSAs (default: `False`) |


## Files output
* [*prefix*.*domain*.filtered.tsv](../files/filtered.tsv.html)
* [*prefix*.*domain*.msa.fasta](../files/msa.fasta.html)
* [*prefix*.*domain*.user_msa.fasta](../files/user_msa.fasta.html)
* [*prefix*.log](../files/gtdbtk.log.html)
* [*prefix*.warnings.log](../files/gtdbtk.warnings.log.html)
* [align/intermediate_results/*prefix*.*domain*.marker_info.tsv](../files/marker_info.tsv.html)

## Example

### Input

```bash
gtdbtk align --identify_dir identify_output/ --out_dir align_output --cpus 3
```

### Output

```text
[2020-04-14 09:14:44] INFO: GTDB-Tk v1.1.0
[2020-04-14 09:14:44] INFO: gtdbtk align --identify_dir identify_output/ --out_dir align_output --cpus 3
[2020-04-14 09:14:44] INFO: Using GTDB-Tk reference data version r89: /release89
[2020-04-14 09:14:44] INFO: Aligning markers in 3 genomes with 3 threads.
[2020-04-14 09:14:44] INFO: Processing 3 genomes identified as archaeal.
[2020-04-14 09:14:44] INFO: Read concatenated alignment for 1248 GTDB genomes.
==> Finished aligning 3 of 3 (100.0%) genomes.
[2020-04-14 09:14:49] INFO: Masking columns of multiple sequence alignment using canonical mask.
[2020-04-14 09:14:52] INFO: Masked alignment from 32675 to 5124 AAs.
[2020-04-14 09:14:52] INFO: 0 user genomes have amino acids in <10.0% of columns in filtered MSA.
[2020-04-14 09:14:52] INFO: Creating concatenated alignment for 1251 GTDB and user genomes.
[2020-04-14 09:14:52] INFO: Creating concatenated alignment for 3 user genomes.
[2020-04-14 09:14:52] INFO: Done.
```