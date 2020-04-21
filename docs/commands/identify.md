---
layout: default
title: identify
nav_order: 4
parent: Commands
---

# identify

Identify marker genes in genome(s).

## Arguments

### Required

| Argument     | Description                |
|:------------|:---------------------------|
| `--out_dir OUT_DIR`          | directory to output files |


#### Mutually exclusive

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
| `--force`          | Try to continue processing if an error occurs on a single genome (default: `False`). |
| `-h`, `--h`           | Display the help message. |


## Files output

* [*prefix*.*domain*.markers_summary.tsv](../files/markers_summary.tsv.html)
* [*prefix*.log](../files/gtdbtk.log.html)
* [*prefix*.translation_table_summary.tsv](../files/translation_table_summary.tsv.html)
* [*prefix*.warnings.log](../files/gtdbtk.warnings.log.html)
* identify/intermediate_results/marker_genes/*genome_id*/
    * [*genome_id*_pfam_tophit.tsv](../files/pfam_tophit.tsv.html)
    * [*genome_id*_pfam.tsv](../files/pfam.tsv.html)
    * [*genome_id*_protein.faa](../files/protein.faa.html)
    * [*genome_id*_protein.fna](../files/protein.fna.html)
    * [*genome_id*_protein.gff](../files/protein.gff.html)
    * [*genome_id*_tigrfam.out](../files/tigrfam.out.html)
    * [*genome_id*_tigrfam_tophit.tsv](../files/tigrfam_tophit.tsv.html)
    * [*genome_id*_tigrfam.tsv](../files/tigrfam.tsv.html)
    * [prodigal_translation_table.tsv](../files/prodigal_translation_table.tsv.html)

## Example

### Input

```bash
gtdbtk identify --genome_dir genomes/ --out_dir identify_output --cpus 3
```

### Output

```text
[2020-04-14 08:51:00] INFO: GTDB-Tk v1.1.0
[2020-04-14 08:51:00] INFO: gtdbtk identify --genome_dir genomes/ --out_dir identify_output --cpus 3
[2020-04-14 08:51:00] INFO: Using GTDB-Tk reference data version r89: /release89
[2020-04-14 08:51:00] INFO: Identifying markers in 3 genomes with 3 threads.
[2020-04-14 08:51:00] INFO: Running Prodigal V2.6.3 to identify genes.
==> Finished processing 3 of 3 (100.0%) genomes.
[2020-04-14 08:51:18] INFO: Identifying TIGRFAM protein families.
==> Finished processing 3 of 3 (100.0%) genomes.
[2020-04-14 08:51:27] INFO: Identifying Pfam protein families.
==> Finished processing 3 of 3 (100.0%) genomes.
[2020-04-14 08:51:29] INFO: Annotations done using HMMER 3.1b2 (February 2015)
[2020-04-14 08:51:29] INFO: Done.
```