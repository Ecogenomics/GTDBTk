---
layout: default
title: ani_rep
nav_order: 10
parent: Commands
---

# ani_rep

Compute the ANI of input genomes to all GTDB-Tk representative genomes.

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

#### Mash arguments

| Argument   | Description                |
|:------------------|:---------------------------|
| `--no_mash`          | Skip pre-filtering using Mash (default: False). |
| `--mash_k MASH_K`          | k-mer size [1-32] (default: `21`). |
| `--mash_s MASH_S`          | maximum number of non-redundant hashes (default: `1000`). |
| `--mash_d MASH_D`          | maximum distance to keep [0-1] (default: `0.3`). |
| `--mash_v MASH_V`          | maximum p-value to keep [0-1]) (default: `1.0`). |

#### FastANI Arguments

| Argument   | Description                |
|:------------------|:---------------------------|
| `--min_af MIN_AF`          | Minimum alignment fraction to consider for genomes (default: `0.65`). |

#### GTDB-Tk Arguments

| Argument   | Description                |
|:------------------|:---------------------------|
| `--extension EXTENSION`          | Extension of files to process, `gz` = gzipped (default: `fna`). |
| `--prefix PREFIX`          | Prefix for output files (default: `gtdbtk`). |
| `--cpus CPUS`          | The number of CPUs to use (default: `1`). |
| `-h`, `--h`           | Display the help message. |


## Files output

* [*prefix*.ani_closest.tsv](../files/ani_closest.tsv.html)
* [*prefix*.ani_summary.tsv](../files/ani_summary.tsv.html)
* [*prefix*.log](../files/gtdbtk.log.html)
* [*prefix*.warnings.log](../files/gtdbtk.warnings.log.html)
* [*prefix*.warnings.log](../files/gtdbtk.warnings.log.html)
* intermediate_results/mash/
    * [*prefix*.gtdb_ref_sketch.msh](../files/gtdbtk_ref_sketch.msh.html)
    * [*prefix*.mash_distances.msh](../files/mash_distances.msh.html)
    * [*prefix*.user_query_sketch.msh](../files/user_query_sketch.msh.html)


## Example

### Input
    
```bash
gtdbtk ani_rep --genome_dir genomes/ --out_dir ani_rep/ --cpus 70
```

### Output

```text
[2020-04-13 10:51:58] INFO: GTDB-Tk v1.1.0
[2020-04-13 10:51:58] INFO: gtdbtk ani_rep --genome_dir genomes/ --out_dir ani_rep/ --cpus 70
[2020-04-13 10:51:58] INFO: Using GTDB-Tk reference data version r89: /release89
[2020-04-13 10:51:59] INFO: Using Mash version 2.2.2
[2020-04-13 10:51:59] INFO: Creating Mash sketch file: ani_rep/intermediate_results/mash/gtdbtk.user_query_sketch.msh
==> Sketching 3 of 3 (100.0%) genomes
[2020-04-13 10:51:59] INFO: Creating Mash sketch file: ani_rep/intermediate_results/mash/gtdbtk.gtdb_ref_sketch.msh
==> Sketching 24706 of 24706 (100.0%) genomes
[2020-04-13 10:53:13] INFO: Calculating Mash distances.
[2020-04-13 10:53:14] INFO: Calculating ANI with FastANI.
==> Processing 874 of 874 (100.0%) comparisons.
[2020-04-13 10:53:23] INFO: Summary of results saved to: ani_rep/gtdbtk.ani_summary.tsv
[2020-04-13 10:53:23] INFO: Closest representative hits saved to: ani_rep/gtdbtk.ani_closest.tsv
[2020-04-13 10:53:23] INFO: Done.
```