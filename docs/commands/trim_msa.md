---
layout: default
title: trim_msa
nav_order: 11
parent: Commands
---

# trim_msa

The `trim_msa` command will trim a MSA given a user-specified mask file, or the archaeal/bacterial 
mask present in the reference data.

## Arguments

### Required

| Argument      | Description                |
|:-----------------|:---------------------------|
| `--untrimmed_msa UNTRIMMED_MSA`        | The path to the untrimmed MSA. |
| `--output OUTPUT`            | The path to save the trimmed MSA. |

#### Mutually exclusive

| Argument     | Description                |
|:------------|:---------------------------|
| `--mask_file MASK_FILE`          | Path to user-specified mask file. |
| `--reference_mask {arc, bac}`        | Use the GTDB-Tk `arc`haeal, or `bac`terial reference mask. |

### Optional

| Argument   | Description                |
|:------------------|:---------------------------|
| `-h`, `--h`           | Display the help message. |


## Example 

### Input

```bash
gtdbtk trim_msa --untrimmed_msa msa.faa --output msa_trim.faa --mask_file mask.txt
```

#### msa.faa
```text
>genome_a
AKLAK
```

#### mask.txt
```text
01011
```

### Output

```
[2020-04-13 10:25:13] INFO: GTDB-Tk v1.1.0
[2020-04-13 10:25:13] INFO: gtdbtk trim_msa --untrimmed_msa msa.faa --output msa_trim.faa --mask_file mask.txt
[2020-04-13 10:25:13] INFO: Using GTDB-Tk reference data version r89: /release89
[2020-04-13 10:25:13] INFO: Done.
```

#### msa_trim.faa

```text
>genome_a
KAK
```

