---
layout: default
title: infer
nav_order: 6
parent: Commands
---

# infer

Infer tree from multiple sequence alignment.

## Arguments

### Required

| Argument     | Description                |
|:------------|:---------------------------|
| `--msa_file MSA_FILE`          | multiple sequence alignment in FASTA format |
| `--out_dir OUT_DIR`          | directory to output files |

### Optional

| Argument   | Description                |
|:------------------|:---------------------------|
| `--prot_model {JTT,WAG,LG}`          | protein substitution model for tree inference (default: `WAG`). |
| `--no_support`          | do not compute local support values using the Shimodaira-Hasegawa test (default: `False`). |
| `--no_gamma`          | do not rescale branch lengths to optimize the Gamma20 likelihood (default: `False`). |
| `--prefix PREFIX`          | desired prefix for output files (default: `gtdbtk`). |
| `--cpus CPUS`          | The number of CPUs to use (default: `1`). |
| `-h`, `--h`           | Display the help message. |

## Files output

* [*prefix*.log](../files/gtdbtk.log.html)
* [*prefix*.unrooted.tree](../files/unrooted.tree.html)
* [*prefix*.warnings.log](../files/gtdbtk.warnings.log.html)
* infer/intermediate_results/
    * [*prefix*.fasttree.log](../files/fasttree.log.html)
    * [*prefix*.tree.log](../files/tree.log.html)


## Example

### Input

```bash
gtdbtk infer --msa_file msa.faa --out_dir infer_out
```

### Output

```text
[2020-04-14 09:37:55] INFO: GTDB-Tk v1.1.0
[2020-04-14 09:37:55] INFO: gtdbtk infer --msa_file msa.faa --out_dir infer_out
[2020-04-14 09:37:55] INFO: Using GTDB-Tk reference data version r89: /release89
[2020-04-14 09:37:55] INFO: Inferring FastTree (WAG, +gamma, support) using a maximum of 1 CPUs.
[2020-04-14 09:37:55] INFO: FastTree version: 2.1.10
[2020-04-14 09:37:55] INFO: FastTree version: 2.1.10
[2020-04-14 09:37:55] INFO: Done.
```