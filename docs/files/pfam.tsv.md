---
layout: default
title: pfam.tsv
nav_order: 14
parent: Files
---

# *genome_id*.pfam.tsv

The raw output of the pfam search.

## Produced by
 * [identify](../commands/identify.html)
 * [classify_wf](../commands/classify_wf.html)
 
## Example
```text
# <seq id> <alignment start> <alignment end> <envelope start> <envelope end> <hmm acc> <hmm name> <type> <hmm start> <hmm end> <hmm length> <bit score> <E-value> <significance> <clan> 
NC_013790.1_1086     29    400     28    400 PF00368.13  HMG-CoA_red       Family     2   373   373    470.1  5.2e-142   1 No_clan  
NC_013790.1_1218     64    245     64    245 PF02006.11  DUF137            Family     1   178   178    232.0   4.3e-70   1 No_clan  
NC_013790.1_1352      2    180      1    182 PF00827.12  Ribosomal_L15e    Family     2   182   192    278.4   3.4e-84   1 No_clan  
NC_013790.1_1468     19    100     19    102 PF01282.14  Ribosomal_S24e    Family     1    82    84     71.1   7.2e-21   1 No_clan  
NC_013790.1_1469     45    165     43    166 PF04019.7   DUF359            Family     3   120   121    117.0   4.3e-35   1 No_clan  
NC_013790.1_1475      3    120      1    121 PF01092.14  Ribosomal_S6e     Family     3   126   127    115.1   2.1e-34   1 No_clan
```