---
layout: default
title: export_msa
nav_order: 12
parent: Commands
---

# export_msa


The `export_msa` will export the untrimmed archaeal or bacterial MSA used in the reference data.

## Arguments

### Required

| Argument       | Description                |
|:----------------|:---------------------------|
| `--domain {arc, bac}`           | The domain to export `arc`haea or `bac`teria. |
| `--output OUTPUT`          | The path to the output file. |

### Optional

| Argument      | Description                |
|:-------------------|:---------------------------|
| `-h`, `--h`           | Display the help message. |

## Example

### Input
```bash
gtdbtk export_msa --domain arc --output /tmp/msa.faa
```

### Output

| Path     |   Description                |
|:-------------|:---------------------------|
| msa.faa    | The untrimmed MSA from the reference data for the specified domain.   |


```
[2020-04-13 10:03:05] INFO: GTDB-Tk v1.1.0
[2020-04-13 10:03:05] INFO: gtdbtk export_msa --domain arc --output /tmp/msa.faa
[2020-04-13 10:03:05] INFO: Using GTDB-Tk reference data version r89: /release89
[2020-04-13 10:03:05] INFO: Done.
```
