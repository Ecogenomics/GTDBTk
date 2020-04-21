---
layout: default
title: root
nav_order: 8
parent: Commands
---

# root

In development
{: .label .label-yellow }

Root a tree using an outgroup.

## Arguments

### Required

| Argument     | Description                |
|:------------|:---------------------------|
| `--input_tree INPUT_TREE`          | tree to root in newick format |
| `--outgroup_taxon OUTGROUP_TAXON`          |taxon to use as outgroup (e.g., `p__Patescibacteria`) |
| `--output_tree OUTPUT_TREE`          |path to write rooted tree |


### Optional

| Argument   | Description                |
|:------------------|:---------------------------|
| `--custom_taxonomy_file`           | custom taxonomy string for at least the genomes belonging to the outgroup |
| `-h`, `--h`           | Display the help message. |


## Example

### Input
    
```bash
gtdbtk root --input_tree input.tree --outgroup_taxon p__Nanoarchaeota --output_tree output.tree
```

### Output

```text
[2020-04-14 08:26:53] INFO: GTDB-Tk v1.1.0
[2020-04-14 08:26:53] INFO: gtdbtk root --input_tree input.tree --outgroup_taxon p__Nanoarchaeota --output_tree output.tree
[2020-04-14 08:26:53] INFO: Using GTDB-Tk reference data version r89: /release89
[2020-04-14 08:26:53] WARNING: Tree rooting is still under development!
[2020-04-14 08:26:53] INFO: Identifying genomes from the specified outgroup.
[2020-04-14 08:26:53] INFO: Identified 101 outgroup taxa in the tree.
[2020-04-14 08:26:53] INFO: Identified 1151 ingroup taxa in the tree.
[2020-04-14 08:26:53] INFO: Outgroup is monophyletic.
[2020-04-14 08:26:53] INFO: Rerooting tree.
[2020-04-14 08:26:53] INFO: Rerooted tree written to: output.tree
[2020-04-14 08:26:53] INFO: Done.
```