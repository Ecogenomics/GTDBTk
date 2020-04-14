---
layout: default
title: decorate
nav_order: 9
parent: Commands
---

# decorate

In development
{: .label .label-yellow }

Decorate a tree with the GTDB-Tk taxonomy.

## Arguments

### Required

| Argument     | Description                |
|:------------|:---------------------------|
| `--input_tree INPUT_TREE`          | tree to decorate in newick format |
| `--output_tree OUTPUT_TREE`          | path to write decorated tree |


### Optional

| Argument   | Description                |
|:------------------|:---------------------------|
| `-h`, `--h`           | Display the help message. |


## Example

### Input
    
```bash
gtdbtk decorate --input_tree input.tree --output_tree output.tree
```

### Output

```text
[2020-04-14 08:20:51] INFO: GTDB-Tk v1.1.0
[2020-04-14 08:20:51] INFO: gtdbtk decorate --input_tree input.tree --output_tree output.tree
[2020-04-14 08:20:51] INFO: Using GTDB-Tk reference data version r89: /release89
[2020-04-14 08:20:51] WARNING: DECORATE NOT YET IMPLEMENTED!
[2020-04-14 08:20:51] INFO: Done.
```