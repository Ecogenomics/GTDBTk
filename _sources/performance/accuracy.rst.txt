.. _performance/Accuracy:


Accuracy
========

The similarity of GTDB-Tk v1 and v2 classifications were first assessed using 16,710 bacterial genomes from the GEMs dataset (Nayfach et al., 2021) that represent novel taxa relative to GTDB R07-RS207.
| Only 12 genomes (0.07%) did not have identical classifications between GTDB-Tk v1 and the divide-and-conquer approach used in GTDB-Tk v2.
| The majority of incongruence was due to genomes being over- (6 genomes) or under-classified (4 genomes) by a single taxonomic rank. Only 2 genomes had conflicting taxonomic assignments, and these were both relatively poor-quality genomes assigned as new classes in alternative phyla.

.. flat-table:: Table 1. Novelty of GEM genomes relative to GTDB R07-RS207 based on GTDB-Tk v1 classifications.
   :header-rows: 2

   * -
     -
     - :cspan:`4` GTDB-Tk v2 classifications relative to GTDB-Tk v1 classifications
   * - Toxon Novelty
     - No genomes
     - Congruent
     - Conflict
     - Underclassified
     - Overclassified
   * - Novel phylum
     - 3
     - 2
     - 0
     - 0
     - 1
   * - Novel class
     - 42
     - 35
     - 2
     - 2
     - 2
   * - Novel order
     - 144
     - 143
     - 0
     - 0
     - 1
   * - Novel family
     - 543
     - 540
     - 0
     - 1
     - 2
   * - Novel genus
     - 3,222
     - 3,219
     - 0
     - 1
     - 0
   * - Novel species
     - 12,756
     - 12,756
     - 0
     - 0
     - 0