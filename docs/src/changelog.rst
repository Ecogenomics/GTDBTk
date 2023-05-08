
Change log
==========

2.3.0
-----

Bug Fixes:

* (`#508 <https://github.com/Ecogenomics/GTDBTk/issues/508>`_) (`#509 <https://github.com/Ecogenomics/GTDBTk/issues/509>`_) If **ALL** genomes for a specific domain are either filtered out or classified with ANI they are now reported in the summary file.

Minor changes:

* (`#491 <https://github.com/Ecogenomics/GTDBTk/issues/491>`_) (`#498 <https://github.com/Ecogenomics/GTDBTk/issues/498>`_) Allow GTDB-Tk to show ``--help`` and ``-v`` without ``GTDBTK_DATA_PATH`` being set.
  * WARNING: This is a breaking change if you are importing GTDB-Tk as a library and importing values from ``gtdbtk.config.config``, instead you need to import as ``from gtdbtk.config.common import CONFIG`` then access values via ``CONFIG.<var>``
* (`#508 <https://github.com/Ecogenomics/GTDBTk/issues/508>`_) Mash distance is changed from 0.1 to 0.15 . This is will increase the number of FastANI comparisons but will cover cases wheere genomes have a larger Mash distance but a small ANI.
* (`#497 <https://github.com/Ecogenomics/GTDBTk/issues/497>`_) Add a ``convert_to_species`` function is GTDB-Tk to replace GCA/GCF ids with their GTDB species name
* Add ``--db_version`` flag to ``check_install`` to check the version of previous GTDB-Tk packages.

2.2.6
-----

Bug Fixes:

* (`#493 <https://github.com/Ecogenomics/GTDBTk/issues/493>`_) Fix issue with --full-tree flag (related to skipping ANI steps)

Minor changes:

* Change URL for documentation to 'https://ecogenomics.github.io/GTDBTk/installing/index.html'
* Improve portability of the ANI_screen step by regenerating the paths of reference genomes in the current filesystem for mash_db.msh


2.2.5
-----

Bug Fixes:

* gtdbtk.json is now reset when the pipeline is re run and the status of ani_Screen is not 'complete'

Minor changes:

* When using '--genes' , ANI steps are skipped and Warnings are raised to the user to
    inform them that classification is less accurate.
* (`#486 <https://github.com/Ecogenomics/GTDBTk/issues/486>`_) Environment variables can be used in GTDBTK_DATA_PATH
* 'is_consistent' function in 'mash.py' compares only the filenames, not the full paths
* Add cutoff arguments to PfamScan ( Thanks @AroneyS for the contribution)


2.2.5
-----

Bug Fixes:

* gtdbtk.json is now reset when the pipeline is re run and the status of ani_Screen is not 'complete'

Minor changes:

* When using '--genes' , ANI steps are skipped and Warnings are raised to the user to inform them that classification is less accurate.
* (`#486 <https://github.com/Ecogenomics/GTDBTk/issues/486>`_) Environment variables can be used in GTDBTK_DATA_PATH
* 'is_consistent' function in 'mash.py' compares only the filenames, not the full paths
* Add cutoff arguments to PfamScan ( Thanks @AroneyS for the contribution)

2.2.4
-----

Bug Fixes:

* (`#475 <https://github.com/Ecogenomics/GTDBTk/issues/475>`_) If all genomes are classified using ANI, Tk will skip the identify step and align steps

Minor changes:

* Add hidden '--skip_pplacer' flag to skip pplacer step ( useful for debugging)
* Improve documentation
* Convert stage_logger to a Singleton class
* Use existing ANI results if available


2.2.3
-----

Bug Fixes:

* Fix prodigal_fail_counter issue

2.2.2
-----

Bug Fixes:

* (`#471 <https://github.com/Ecogenomics/GTDBTk/issues/471>`_) Fix Pplacer issue


2.2.1
-----

Build:

* (`#470 <https://github.com/Ecogenomics/GTDBTk/issues/470>`_) Add missing Pydantic dependency.


2.2.0
-----

Minor changes:

* (`#433 <https://github.com/Ecogenomics/GTDBTk/issues/433>`_) Added additional checks to ensure that the `--outgroup_taxon` cannot be set to a domain (`root`, `de_novo_wf`).
* (`#459 <https://github.com/Ecogenomics/GTDBTk/issues/459>`_ / `#462 <https://github.com/Ecogenomics/GTDBTk/issues/462>`_ ) Fix deprecated np.bool in prodigal_biolib.py. Special thanks to @neoformit for his contribution.
* (`#466 <http://github.com/Ecogenomics/GTDBTk/issues/466>`_) RED value has been rounded to 5 decimals after the comma.
* (`#451 <http://github.com/Ecogenomics/GTDBTk/issues/451>`_) Extra checks have been added when Prodigal fails.
* (`#448 <http://github.com/Ecogenomics/GTDBTk/issues/448>`_) Warning has been added when all the genomes are filtered out and not classified.

Bug Fixes:

* (`#420 <https://github.com/Ecogenomics/GTDBTk/issues/420>`_) Fixed an issue where GTDB-Tk might hang when classifying TIGRFAM markers (`identify`, `classify_wf`, `de_novo_wf`). Special thanks to @lfenske-93 and @sjaenick for their contribution.
* (`#428 <https://github.com/Ecogenomics/GTDBTk/issues/428>`_) Fixed an issue where the `--gtdbtk_classification_file` would raise an error trying to read the `classify` summary (`root`, `de_novo_wf`).
* (`#439 <https://github.com/Ecogenomics/GTDBTk/issues/439>`_) Fix the pipeline when using protein files instead of nucleotide files. symlink uses absolute path instead.




2.1.1
-----

Documentation:

* (`#410 <https://github.com/Ecogenomics/GTDBTk/issues/410>`_) Add documentation for `convert_to_itol`

Bug Fixes:

* (`#399 <https://github.com/Ecogenomics/GTDBTk/issues/399>`_) Fix `--genes` option attempting to create a directory.
* (`#400 <https://github.com/Ecogenomics/GTDBTk/issues/400>`_) Updated contig.py to fix inconsistent pplacer paths causing the program to crash.


2.1.0
-----

Major changes:

* GTDB-TK now uses a **divide-and-conquer** approach where the bacterial reference tree is split into multiple **class**-level subtrees. This reduces the memory requirements of GTDB-Tk from **320 GB** of RAM when using the full GTDB R07-RS207 reference tree to approximately **55 GB**. A manuscript describing this approach is in preparation. If you wish to continue using the full GTDB reference tree use the `--full-tree` flag. This is the main change from v2.0.0. The split tree approach has been modified from order-level trees to class-level trees to resolve specific classification issues (see `#383 <https://github.com/Ecogenomics/GTDBTk/issues/383>`_).
* Genomes that cannot be assigned to a domain (e.g. genomes with no bacterial or archaeal markers or genomes with no genes called by Prodigal) are now reported in the `gtdbtk.bac120.summary.tsv` as 'Unclassified'
* Genomes filtered out during the alignment step are now reported in the `gtdbtk.bac120.summary.tsv` or `gtdbtk.ar53.summary.tsv` as 'Unclassified Bacteria/Archaea'
* `--write_single_copy_genes` flag in now available in the `classify_wf` and `de_novo_wf` workflows.


Features:

* (`#392 <https://github.com/Ecogenomics/GTDBTk/issues/392>`_) `--write_single_copy_genes` flag available in workflows.
* (`#387 <https://github.com/Ecogenomics/GTDBTk/issues/392>`_) specific memory requirements set in classify_wf depending on the classification approach.


2.0.0
-----

Major changes:

* GTDB-TK now uses a **divide-and-conquer** approach where the bacterial reference tree is split into multiple order-level subtrees. This reduces the memory requirements of GTDB-Tk from **320 GB** of RAM when using the full GTDB R07-RS207 reference tree to approximately **35 GB**. A manuscript describing this approach is in preparation. If you wish to continue using the full GTDB reference tree use the `--full-tree` flag.
* Archaeal classification now uses a refined set of 53 archaeal-specific marker genes based on the recent publication by `Dombrowski et al., 2020 <https://www.nature.com/articles/s41467-020-17408-w>`_. This set of archaeal marker genes is now used by GTDB for curating the archaeal taxonomy.
* By default, all directories containing intermediate results are **now removed** by default at the end of the `classify_wf` and `de_novo_wf` pipelines. If you wish to retain these intermediates files use the `--keep-intermediates` flag.
* All MSA files produced by the `align` step are now compressed with gzip.
* The classification summary and failed genomes files are now the only files linked in the root directory of `classify_wf`.


Features:

* (`#373 <https://github.com/Ecogenomics/GTDBTk/issues/373>`_) `convert_to_itol` to convert trees into iTOL format
* (`#369 <https://github.com/Ecogenomics/GTDBTk/issues/369>`_) Output FASTA files are compressed by default
* (`#369 <https://github.com/Ecogenomics/GTDBTk/issues/369>`_) Intermediate files will be removed by default when using classify/de-novo workflows unless specified by `--keep_intermediates`
* (`#362 <https://github.com/Ecogenomics/GTDBTk/issues/362>`_) Add --genes flag for Error
* (`#360 <https://github.com/Ecogenomics/GTDBTk/issues/360>`_ / `#356 <https://github.com/Ecogenomics/GTDBTk/issues/356>`_) A warning will be displayed if pplacer fails to place a genome

**Important**

* This version is **not** backwards compatible with GTDB release 202.
* This version requires a `new reference package <https://data.gtdb.ecogenomic.org/releases/release207/207.0/auxillary_files/gtdbtk_r207_data.tar.gz>`_


1.7.0
-----

* (`#336 <https://github.com/Ecogenomics/GTDBTk/issues/336>`_) Warn the user if they have provided an incorrectly formatted taxonomy file.
* (`#348 <https://github.com/Ecogenomics/GTDBTk/issues/348>`_) Gracefully exit the program if no single copy hits could be identified.
* (`#351 <https://github.com/Ecogenomics/GTDBTk/issues/351>`_) Fixed an issue where GTDB-Tk would crash if spaces were present in the reference data path.
* (`#354 <https://github.com/Ecogenomics/GTDBTk/pull/354>`_) Added optional ``--tmpdir`` argument to set temporary directory (thanks `tr11-sanger <https://github.com/tr11-sanger>`_!).


1.6.0
-----

* (`#337 <https://github.com/Ecogenomics/GTDBTk/issues/337>`_) Set minimum `tqdm` version to `4.35.0`
* (`#335 <https://github.com/Ecogenomics/GTDBTk/pull/335>`_) Fixed typo in output log messages (@fplaza)
* Removed the option to re-calculate RED values (`--recalculate_red`)

1.5.1
-----

* (`#327 <https://github.com/Ecogenomics/GTDBTk/issues/327>`_) Disallow spaces in genome names/file paths due to downstream application issues.
* (`#326 <https://github.com/Ecogenomics/GTDBTk/issues/326>`_) Disallow genome names that are blank.

1.5.0
-----

* (`#311 <https://github.com/Ecogenomics/GTDBTk/issues/311>`_) Updated GTDB-Tk to support R202.
  See https://ecogenomics.github.io/GTDBTk/installing/index.html#gtdb-tk-reference-data for instructions on downloading R202.


1.4.2
-----

* (`#311 <https://github.com/Ecogenomics/GTDBTk/issues/311>`_) Fixed --scratch_dir not working in v 1.4.1 for classify_wf
* (`#312 <https://github.com/Ecogenomics/GTDBTk/issues/311>`_) Automatic drop of genome leads to error in downstream modules of classify_wf


1.4.1
-----

* Updated GitHub CI/CD to trigger docker build / tag version on release.
* (`#255 <https://github.com/Ecogenomics/GTDBTk/issues/255>`_) (`#297 <https://github.com/Ecogenomics/GTDBTk/issues/297>`_)
  Fixed ``'Namespace' object has no attribute`` errors by adding default arguments to argparse.


1.4.0
-----

* Check if stdout is being piped to a file before adding colour.
* (`#283 <https://github.com/Ecogenomics/GTDBTk/issues/283>`_) Significantly improved ``classify`` performance (noticeable when running trees > 1,000 taxa).
* Automatically cap pplacer CPUs to 64 unless specifying ``--pplacer_cpus`` to prevent pplacer from hanging.
* (`#262 <https://github.com/Ecogenomics/GTDBTk/issues/262>`_) Added ``--write_single_copy_genes`` to the ``identify`` command. Writes unaligned single-copy AR53/BAC120 marker genes to disk.
* When running ``-version`` warn if GTDB-Tk is not running the most up-to-date version (disable via ``GTDBTK_VER_CHECK = False`` in ``config.py``). If GTDB-Tk encounters an error it will silently continue (3 second timeout).
* (`#276 <https://github.com/Ecogenomics/GTDBTk/issues/276>`_) Renamed the column ``aa_percent`` to ``msa_percent`` in ``summary.tsv`` (produced by ``classify``).
* (`#286 <https://github.com/Ecogenomics/GTDBTk/pull/286>`_) Fixed a file not found error when the reference data is a symbolic link (thanks `davidealbanese <https://github.com/davidealbanese>`_!).
* (`#277 <https://github.com/Ecogenomics/GTDBTk/issues/277>`_) Fixed an issue where if the user overrides the translation table using the optional 3rd column in the batchfile, the other coding density would appear as -100. Both translation table densities are now reported.
* The :ref:`commands/check_install` command now also checks that all third party binaries can be found on the system path.
* The ``align`` step is now approximately 10x faster.
* (`#289 <https://github.com/Ecogenomics/GTDBTk/issues/289>`_) Added ``--min_af`` to ``classify`` and ``classify_wf`` which allows the user to specify the minimum alignment fraction for FastANI.
* Added the ``--mash_db`` command to re-use the GTDB-Tk Mash reference database in ``ani_rep``.


1.3.0
-----

* This version of GTDB-Tk requires a new version of the GTDB-Tk reference package
  (gtdbtk_r95_data.tar.gz) `available here <https://data.ace.uq.edu.au/public/gtdb/data/releases/release95/95.0/auxillary_files/gtdbtk_r95_data.tar.gz>`_.
* Updated reference package to use the **GTDB Release 95** taxonomy.
* Report if the species-specific ANI circumscription criteria is satisfied in the ``ani_closest.tsv`` file output by ``ani_rep``.
* Estimated time until completion has been dampened.


1.2.0
-----

* (`#241 <https://github.com/Ecogenomics/GTDBTk/issues/241>`_) Moved GTDB-Tk entry point to ``__main__.py`` instead of
  ``bin/gtdbtk`` to support execution in some HPC systems (``gtdbtk`` will still be aliased on install).
* (`#251 <https://github.com/Ecogenomics/GTDBTk/issues/251>`_) Allow parsing of FastANI v1.0 output files. However, a warning will be displayed to update FastANI.
* (`#254 <https://github.com/Ecogenomics/GTDBTk/issues/254>`_) Fixed an issue where ``--scratch_dir`` would fail, and not clean-up the mmap file.
* (`#242 <https://github.com/Ecogenomics/GTDBTk/pull/242>`_) Added the ``decorate`` command allowing the ``de novo workflow`` to be run
* (`#244 <https://github.com/Ecogenomics/GTDBTk/pull/244>`_) Added the ``infer_rank`` method which established the taxonomic ranks of internal nodes of user trees based on RED
* (`#248 <https://github.com/Ecogenomics/GTDBTk/pull/248>`_) If the identify command is run on the same directory, genomes which were already processed will be skipped.
* (`#248 <https://github.com/Ecogenomics/GTDBTk/pull/248>`_) Improved ``pplacer`` output with running the ``classify`` command.


1.1.0
-----

* In rare cases pplacer would assign an empty taxonomy string which would raise an error.
* (`#229 <https://github.com/Ecogenomics/GTDBTk/issues/229>`_) Genomes using windows line carriage ``\r\n`` would raise an error.
* (`#227 <https://github.com/Ecogenomics/GTDBTk/issues/227>`_) CentOS machines would fail when using ``~`` in paths.
* The bac120 symlink was pointing to the archaeal tree when using the ``root`` command.
* Updated the ``gtdb_to_ncbi_majority_vote.py`` script for translating taxonomy.
* (`#195 <https://github.com/Ecogenomics/GTDBTk/issues/195>`_) Added the ``--pplacer_cpus`` argument to specify the number of pplacer threads when running ``classify`` and ``classify_wf`` (#195).
* (`#198 <https://github.com/Ecogenomics/GTDBTk/issues/198>`_) The ``--debug`` flag of ``align`` outputs aligned markers to disk before trimming.
* (`#225 <https://github.com/Ecogenomics/GTDBTk/issues/225>`_) An optional third column in the ``--batchfile`` will specify an override to which translation table should be used.
  Leave blank to automatically determine the translation table (default).
* (`#131 <https://github.com/Ecogenomics/GTDBTk/issues/131>`_) Users can now specify genomes which have NCBI accessions, as long as they are not GTDB-Tk
  representatives (a warning will be raised).
* (`#191 <https://github.com/Ecogenomics/GTDBTk/issues/191>`_) Added a new command ``ani_rep`` which calculates the ANI of input genomes to all GTDB
  representative genomes.
* This command uses `Mash <https://github.com/marbl/Mash>`_ in a pre-filtering step. If pre-filtering is enabled (default)
  then ``mash`` will need to be on the system path. To disable pre-filtering use the ``--no_mash`` flag.
* (`#230 <https://github.com/Ecogenomics/GTDBTk/issues/235>`_) Improved how markers are used in determining the correct domain, and gene selection for the alignment.


1.0.2
-----

* Fixed an issue where FastANI threads would timeout with ``FastANI returned a non-zero exit code.``
* Versions affected: ``1.0.0``, and ``1.0.1``.


1.0.0
-----

* Migrated to **Python 3**, you must be running at least **Python 3.6** or later to use this version.
* ``check_install`` now does an exhaustive check of the reference data.
* Resolved an issue where gene calling would fail for low quality genomes (#192).
* Improved FastANI multiprocessing performance.
* Third party software versions are reported where possible.


0.3.3
-----

* A bug has been fixed which affected ``classify`` and ``classify_wf`` when using the ``--batchfile``
  argument with genome IDs that differed from the FASTA filename. This issue resulted in
  the assigned taxonomy being derived only from tree placement without any ANI
  calculations being considered. Consequently, in some cases genomes may have been classified as a new
  species within a genus when they should have been assigned to an existing species. If you have genomes
  with species assignments this bug did not impact you.
* Progress is now displayed for: hmmalign, and pplacer.
* Fixed an issue where the ``root`` command could not be run independently.
* Improved MSA masking performance.


0.3.2
-----

* FastANI calculations are more robust.
* Optimisation of RED calculations.
* Improved output messages when errors are encountered.


0.3.1
-----

* Pplacer taxonomy is now available in the summary file.
* FastANI species assignment will be selected over phylogenetic placement (Topology case).


0.3.0
-----

* Best translation table displayed in summary file.
* GTDB-Tk now supports gzipped genomes as inputs (``--extension gz``).
* By default, GTDB-Tk uses precalculated RED values.
* New option to recalculate RED value during classify step (``--recalculate_red``).
* New option to export the untrimmed reference MSA files.
* New option to skip_trimming during align step.
* New option to use a custom taxonomy file when rooting a tree.
* New FAQ page available.
* New output structure.


0.2.1
-----

* Species classification is now based strictly on the ANI to reference genomes
* The "classify" function now reports the closest reference genome in the summary file even if the ANI is <95%
* The summary.tsv file has 4 new columns: aa_percent, red_values, fastani_reference_radius, and warnings
* By default, the "align" function now performs the same MSA trimming used by the GTDB
* New pplacer support for writing to a scratch file (``--mmap-file`` option)
* Random seed option for MSA trimming has been added to allow for reproducible results
* Configuration of the data directory is now set using the environment variable ``GTDBTK_DATA_PATH`` (see pip installation)
* Perl dependencies has been removed
* Python libraries biolib, mpld3 and jinja have been removed
* This version requires a new version of the GTDB-Tk data package (gtdbtk.r86_v2_data.tar.gz) available `here <https://data.ace.uq.edu.au/public/gtdbtk/release_86/>`_


0.1.3
-----

* GTDB-Tk v0.1.3 has been released and addresses an issue with species assignments based on the placement of
  genomes in the reference tree. This impacted species assignment when submitting multiple closely related genomes.
  Species assignments reported by ANI were not impacted.


0.1.0
-----

* Updated to R86, requires `release 86 <https://data.ace.uq.edu.au/public/gtdbtk/release_86/>`_ to run.
