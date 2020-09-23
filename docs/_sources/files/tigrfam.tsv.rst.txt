.. _files/tigrfam.tsv:

tigrfam.tsv
===========

The raw output produced by Tigrfam parsed into tsv.

Produced by
-----------
 * :ref:`commands/identify`
 * :ref:`commands/classify_wf`


Example
-------

.. code-block:: text

    #                                                               --- full sequence ---- --- best 1 domain ---- --- domain number estimation ----
    # target name        accession  query name           accession    E-value  score  bias   E-value  score  bias   exp reg clu  ov env dom rep inc description of target
    #------------------- ---------- -------------------- ---------- --------- ------ ----- --------- ------ -----   --- --- --- --- --- --- --- --- ---------------------
    NC_013790.1_950      -          TIGR00021            TIGR00021      3e-82  272.2   0.5   3.4e-82  272.0   0.5   1.0   1   0   0   1   1   1   1 # 1219054 # 1219722 # 1 # ID=1_950;partial=00;start_type=ATG;rbs_motif=None;rbs_spacer=None;gc_cont=0.383
    NC_013790.1_1725     -          TIGR00037            TIGR00037    8.2e-50  165.0   0.7   9.2e-50  164.8   0.7   1.0   1   0   0   1   1   1   1 # 2252948 # 2253355 # -1 # ID=1_1725;partial=00;start_type=ATG;rbs_motif=GGAG/GAGG;rbs_spacer=5-10bp;gc_cont=0.348
    NC_013790.1_2080     -          TIGR00042            TIGR00042    4.4e-67  222.0   0.0   4.9e-67  221.8   0.0   1.0   1   0   0   1   1   1   1 # 2740928 # 2741485 # -1 # ID=1_2080;partial=00;start_type=ATG;rbs_motif=None;rbs_spacer=None;gc_cont=0.360
    NC_013790.1_400      -          TIGR00064            TIGR00064   7.6e-102  337.2   5.8  7.6e-102  337.2   5.8   3.4   2   1   2   4   4   1   1 # 571391 # 573106 # 1 # ID=1_400;partial=00;start_type=TTG;rbs_motif=GGA/GAG/AGG;rbs_spacer=5-10bp;gc_cont=0.382

