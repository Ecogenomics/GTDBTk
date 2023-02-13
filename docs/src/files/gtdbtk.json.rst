.. _files/gtdbtk.json:

gtdbtk.json
==========

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
        "version": "2.1.1",
        "command_line": "gtdbtk classify_wf --batchfile /srv/projects/gtdbtk/test_new_features/gems_benchmark/3lines_batchfile.tsv --out_dir /srv/projects/gtdbtk/test_new_features/gems_benchmark/classify_wf_outdir_prescreen_3lines/ --keep_intermediates --cpus 20 --mash_db /srv/projects/gtdbtk/test_new_features/gems_benchmark/mash_sketch/cli/",
        "database_version": "r207",
        "database_path": "/srv/projects/gtdbtk/test_new_features/release207_v2/",
        "steps": [
            {
                "name": "ANI screen",
                "input": "/srv/projects/gtdbtk/test_new_features/gems_benchmark/3lines_batchfile.tsv",
                "output_dir": "/srv/projects/gtdbtk/test_new_features/gems_benchmark/classify_wf_outdir_prescreen_3lines/",
                "output_files": {
                    "bac120": "/srv/projects/gtdbtk/test_new_features/gems_benchmark/classify_wf_outdir_prescreen_3lines/classify/ani_screen/gtdbtk.bac120.ani_summary.tsv"
                },
                "starts_at": "2023-02-01T08:02:17.814231",
                "ends_at": "2023-02-01T08:02:27.782442",
                "duration": "0:00:09",
                "status": "completed",
                "mash_k": 16,
                "mash_s": 5000,
                "mash_v": 1.0,
                "mash_max_dist": 0.1,
                "mash_db": "/srv/projects/gtdbtk/test_new_features/gems_benchmark/mash_sketch/cli/"
            },
