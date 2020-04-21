---
layout: default
title: check_install
nav_order: 14
parent: Commands
---

# check_install

The `check_install` command is used to verify the integrity of the GTDB-Tk reference data.

If any inconsistencies are identified then the program will exit with code 1 and 
`HASH MISMATCH` will be displayed next to the inconsistent item.

## Example
### Input
```bash 
gtdbtk check_install
```

### Output
```
[2020-04-13 09:35:16] INFO: GTDB-Tk v1.1.0
[2020-04-13 09:35:16] INFO: gtdbtk check_install
[2020-04-13 09:35:16] INFO: Using GTDB-Tk reference data version r89: /release89
[2020-04-13 09:35:16] INFO: Running install verification
[2020-04-13 09:35:16] INFO: Checking /release89
[2020-04-13 09:35:16] INFO:          |-- pplacer          OK
[2020-04-13 09:35:16] INFO:          |-- masks            OK
[2020-04-13 09:35:17] INFO:          |-- markers          OK
[2020-04-13 09:35:17] INFO:          |-- radii            OK
[2020-04-13 09:35:20] INFO:          |-- msa              OK
[2020-04-13 09:35:20] INFO:          |-- metadata         OK
[2020-04-13 09:35:20] INFO:          |-- taxonomy         OK
[2020-04-13 09:47:36] INFO:          |-- fastani          OK
[2020-04-13 09:47:36] INFO:          |-- mrca_red         OK
[2020-04-13 09:47:36] INFO: Done.
```