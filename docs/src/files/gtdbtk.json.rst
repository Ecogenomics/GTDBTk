.. _files/gtdbtk.json:

gtdbtk.json
===========

The console output of GTDB-Tk saved to disk in a JSON format.

Produced by
-----------

* :ref:`commands/align`
* :ref:`commands/align`
* :ref:`commands/classify`
* :ref:`commands/classify_wf`
* :ref:`commands/de_novo_wf`
* :ref:`commands/identify`
* :ref:`commands/infer`

Example
-------

.. code-block:: text

    {
    "version": "2.5.0",
    "command_line": "gtdbtk classify_wf --batchfile genomes/3_batchfile.tsv --out_dir classify_wf_3_genomes --cpus 20",
    "database_version": "r226",
    "database_path": "/srv/db/gtdbtk/official/release226",
    "steps": [
        {
            "name": "ANI screen",
            "output_dir": "classify_wf_3_genomes",
            "starts_at": "2025-08-05T20:31:02.156781",
            "ends_at": "2025-08-05T20:35:27.198873",
            "duration": "0:04:25",
            "status": "completed",
            "genome_dir": null,
            "batchfile": "genomes/3_batchfile.tsv",
            "min_af": 0.5,
            "output_files": {}
        },
        {
            "name": "identify",
            "output_dir": "classify_wf_3_genomes",
            "starts_at": "2025-08-05T20:35:27.270095",
            "ends_at": "2025-08-05T20:36:15.865020",
            "duration": "0:00:48",
            "status": "completed",
            "genes": false,
            "extension": "fna",
            "write_single_copy_genes": false,
            "genome_dir": null,
            "batchfile": "genomes/3_batchfile.tsv",
            "output_files": {
                "all": [
                    "classify_wf_3_genomes/identify/gtdbtk.failed_genomes.tsv",
                    "classify_wf_3_genomes/identify/gtdbtk.translation_table_summary.tsv"
                ],
                "ar53": [
                    "classify_wf_3_genomes/identify/gtdbtk.ar53.markers_summary.tsv"
                ],
                "bac120": [
                    "classify_wf_3_genomes/identify/gtdbtk.bac120.markers_summary.tsv"
                ]
            }
        },
        {
            "name": "align",
            "output_dir": "classfiy_wf_3_genomes",
            "starts_at": "2025-08-05T20:36:15.865242",
            "ends_at": "2025-08-05T20:41:38.376388",
            "duration": "0:05:22",
            "status": "completed",
            "identify_dir": "classfiy_wf_3_genomes",
            "skip_gtdb_refs": false,
            "taxa_filter": null,
            "min_perc_aa": 10,
            "custom_msa_filters": false,
            "skip_trimming": false,
            "rnd_seed": null,
            "cols_per_gene": null,
            "min_consensus": null,
            "max_consensus": null,
            "min_perc_taxa": null,
            "outgroup_taxon": null,
            "output_files": {
                "bac120": [
                    "classfiy_wf_3_genomes/align/gtdbtk.bac120.filtered.tsv",
                    "classfiy_wf_3_genomes/align/gtdbtk.bac120.msa.fasta",
                    "classfiy_wf_3_genomes/align/gtdbtk.bac120.user_msa.fasta"
                ]
            }
        },
        {
            "name": "classify",
            "output_dir": "classfiy_wf_3_genomes",
            "starts_at": "2025-08-05T20:41:38.377084",
            "ends_at": "2025-08-05T21:08:25.704484",
            "duration": "0:26:47",
            "status": "completed",
            "align_dir": "classfiy_wf_3_genomes",
            "genome_dir": null,
            "batchfile": "genomes/3_batchfile.tsv",
            "scratch_dir": null,
            "debug_option": false,
            "full_tree": false,
            "skip_ani_screen": true,
            "output_files": {
                "bac120": [
                    "classify_wf_3_genomes/classify/gtdbtk.backbone.bac120.classify.tree",
                    "classify_wf_3_genomes/classify/gtdbtk.bac120.classify.tree.7.tree",
                    "classify_wf_3_genomes/classify/gtdbtk.bac120.classify.tree.2.tree",
                    "classify_wf_3_genomes/classify/gtdbtk.bac120.tree.mapping.tsv",
                    "classify_wf_3_genomes/classify/gtdbtk.bac120.summary.tsv"
                ]
            }
        }
    ],
    "output_dir": "classify_wf_3_genomes",
    "path": "classify_wf_3_genomes/gtdbtk.json"
