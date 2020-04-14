---
layout: default
title: tigrfam.out
nav_order: 25
parent: Files
---

# *genome_id*_tigrfam.out

The raw output produced by Tigrfam.

## Produced by
 * [identify](../commands/identify.html)
 * [classify_wf](../commands/classify_wf.html)


## Example

```text
Query:       TIGR00006  [M=310]
Accession:   TIGR00006
Description: TIGR00006: 16S rRNA (cytosine(1402)-N(4))-methyltransferase
Scores for complete sequences (score includes all domains):
   --- full sequence ---   --- best 1 domain ---    -#dom-
    E-value  score  bias    E-value  score  bias    exp  N  Sequence Description
    ------- ------ -----    ------- ------ -----   ---- --  -------- -----------

   [No hits detected that satisfy reporting thresholds]


Domain annotation for each sequence:

   [No targets detected that satisfy reporting thresholds]


Internal pipeline statistics summary:
-------------------------------------
Query model(s):                              1  (310 nodes)
Target sequences:                         2199  (764686 residues searched)
Passed MSV filter:                       149  (0.0677581); expected 44.0 (0.02)
Passed bias filter:                       65  (0.0295589); expected 44.0 (0.02)
Passed Vit filter:                         5  (0.00227376); expected 2.2 (0.001)
Passed Fwd filter:                         1  (0.000454752); expected 0.0 (1e-05)
Initial search space (Z):               2199  [actual number of targets]
Domain search space  (domZ):               0  [number of targets reported over threshold]
# CPU time: 0.04u 0.00s 00:00:00.04 Elapsed: 00:00:00.04
# Mc/sec: 5926.32
//
Query:       TIGR00019  [M=361]
Accession:   TIGR00019
Description: prfA: peptide chain release factor 1
Scores for complete sequences (score includes all domains):
   --- full sequence ---   --- best 1 domain ---    -#dom-
    E-value  score  bias    E-value  score  bias    exp  N  Sequence Description
    ------- ------ -----    ------- ------ -----   ---- --  -------- -----------

   [No hits detected that satisfy reporting thresholds]


Domain annotation for each sequence:

   [No targets detected that satisfy reporting thresholds]
```