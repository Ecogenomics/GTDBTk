---
layout: default
title: classify
nav_order: 7
parent: Commands
---

# classify

Determine taxonomic classification of genomes.

## Arguments

### Required

| Argument     | Description                |
|:------------|:---------------------------|
| `--align_dir ALIGN_DIR`          | output directory of `align` command |
| `--out_dir OUT_DIR`          | directory to output files |

####  Mutually exclusive 

| Argument     | Description                |
|:------------|:---------------------------|
| `--genome_dir GENOME_DIR`          | directory containing genome files in FASTA format |
| `--batchfile BATCHFILE`        | file describing genomes - tab separated in 2 columns (FASTA file, genome ID) |

### Optional

| Argument   | Description                |
|:------------------|:---------------------------|
| `--extension EXTENSION`          | Extension of files to process, `gz` = gzipped (default: `fna`). |
| `--prefix PREFIX`          | Prefix for output files (default: `gtdbtk`). |
| `--cpus CPUS`          | The number of CPUs to use (default: `1`). |
| `--pplacer_cpus PPLACER_CPUS`          | The number of CPUs to use with pplacer (default: `CPUS`). |
| `--scratch_dir SCRATCH_DIR`          | reduce memory usage by writing to disk (slower). |
| `-r`, `--recalculate_red`          | recalculate RED values based on the reference tree and all added user genomes (default: `False`). |
| `--debug`           | create intermediate files for debugging purposes (default: `False`). |
| `-h`, `--h`           | Display the help message. |



## Files output

* classify
    * intermediate_results
        * [*prefix*.*domain*.classification_pplacer.tsv](../files/classification_pplacer.tsv.html)
        * [*prefix*.*domain*.red_dictionary.tsv](../files/red_dictionary.tsv.html)
        * pplacer
            * [pplacer.*domain*.json](../files/pplacer.domain.json.html)
            * [pplacer.*domain*.out](../files/pplacer.domain.out.html)
* [*prefix*.*domain*.classify.tree](../files/classify.tree.html)
* [*prefix*.*domain*.summary.tsv](../files/summary.tsv.html)
* [*prefix*.log](../files/gtdbtk.log.html)
* [*prefix*.warnings.log](../files/gtdbtk.warnings.log.html)

## Example

### Input

```bash
 gtdbtk classify --genome_dir genomes/ --align_dir align_output/ --out_dir classify_output --cpus 3
```

### Output

```text
[2020-04-14 09:22:35] INFO: GTDB-Tk v1.1.0
[2020-04-14 09:22:35] INFO: gtdbtk classify --genome_dir genomes/ --align_dir align_output/ --out_dir classify_output --cpus 3
[2020-04-14 09:22:35] INFO: Using GTDB-Tk reference data version r89: /release89
[2020-04-14 09:22:35] INFO: Placing 3 archaeal genomes into reference tree with pplacer using 3 cpus (be patient).
Placing genomes |##################################################| 3/3 (100.00%)
[2020-04-14 09:23:32] INFO: pplacer version: v1.1.alpha19-0-g807f6f3
[2020-04-14 09:23:33] INFO: Calculating average nucleotide identity using FastANI.
[2020-04-14 09:23:33] INFO: fastANI version: 1.3
==> Processing 24 of 24 (100.0%) comparisons.
[2020-04-14 09:23:38] INFO: 3 genome(s) have been classified using FastANI and pplacer.
[2020-04-14 09:23:38] INFO: Calculating RED values based on reference tree.
[2020-04-14 09:23:38] INFO: Done.
```