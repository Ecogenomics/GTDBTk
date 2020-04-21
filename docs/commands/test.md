---
layout: default
title: test
nav_order: 13
parent: Commands
---

# test

The `test` command is used to run three small archaeal genomes through the classify workflow.

If your installation is unable to run the `test` command with an exit code of 0, then 
there is an issue with your installation.

## Arguments

### Required

| Argument     | Description                |
|:-----------------|:---------------------------|
| `--out_dir OUT_DIR`      | Directory to output files. |
| `--cpus CPUS`          | The number of CPUs to use. |


### Optional

| Argument      | Description                |
|:-----------------|:---------------------------|
| `-h`, `--h`           | Display the help message. |

## Files output

* [*prefix*.warnings.log](../files/gtdbtk.warnings.log.html)
* [*prefix*.warnings.log](../files/gtdbtk.warnings.log.html)
* [output/](../commands/classify_wf.html#files-output)
* [test_execution.log](../files/test_execution.log.html)

## Example

### Input
```bash
gtdbtk test --out_dir /tmp/test --cpus 3
```

### Output

```
[2020-04-13 09:50:58] INFO: GTDB-Tk v1.1.0
[2020-04-13 09:50:58] INFO: gtdbtk test --out_dir /tmp/test --cpus 3
[2020-04-13 09:50:58] INFO: Using GTDB-Tk reference data version r89: /release89
[2020-04-13 09:50:58] INFO: Command: gtdbtk classify_wf --genome_dir /tmp/test/genomes --out_dir /tmp/test/output --cpus 3
[2020-04-13 09:52:35] INFO: Test has successfully finished.
```