
Change log
==========

1.4.0
-----

* Check if stdout is being piped to a file before adding colour.
* (`#283 <https://github.com/Ecogenomics/GTDBTk/issues/283>`_) Significantly improved ``classify`` performance (noticeable when running trees > 1,000 taxa).
* Automatically cap pplacer CPUs to 64 unless specifying ``--pplacer_cpus`` to prevent pplacer from hanging.
* (`#262 <https://github.com/Ecogenomics/GTDBTk/issues/262>`_) Added ``--write_single_copy_genes`` to the ``identify`` command. Writes unaligned single-copy AR122/BAC120 marker genes to disk.
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
