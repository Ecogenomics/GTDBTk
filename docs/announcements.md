---
layout: default
title: Announcements
nav_order: 1
---


# Announcements

**Note (x, 2020)**:
 * GTDB-Tk v1.1.1 has been released 
   * *Bug fixes:*
        * none?
    * *Features:*
        * Added the `infer_rank` method which established the taxonomic ranks of internal nodes of user trees based on RED
        * Added the `decorate` command allowing the `de novo workflow` to be run

**Note (Apr 9, 2020)**:
* GTDB-Tk v1.1.0 has been released (**we recommend all users update to this version**)
    * *Bug fixes:*
        * In rare cases pplacer would assign an empty taxonomy string which would raise an error.
        * ([#229](https://github.com/Ecogenomics/GTDBTk/issues/229)) Genomes using windows line carriage `\r\n` would raise an error.
        * ([#227](https://github.com/Ecogenomics/GTDBTk/issues/227)) CentOS machines would fail when using `~` in paths.
        * The bac120 symlink was pointing to the archaeal tree when using the `root` command.
    * *Features:*
        * Updated the `gtdb_to_ncbi_majority_vote.py` script for translating taxonomy.
        * ([#195](https://github.com/Ecogenomics/GTDBTk/issues/195)) Added the `--pplacer_cpus` argument to specify the number of pplacer threads when running `classify` and `classify_wf` (#195).
        * ([#198](https://github.com/Ecogenomics/GTDBTk/issues/198)) The `--debug` flag of `align` outputs aligned markers to disk before trimming.
        * ([#225](https://github.com/Ecogenomics/GTDBTk/issues/225)) An optional third column in the `--batchfile` will specify an override to which translation table should be used. 
          Leave blank to automatically determine the translation table (default).
        * ([#131](https://github.com/Ecogenomics/GTDBTk/issues/131)) Users can now specify genomes which have NCBI accessions, as long as they are not GTDB-Tk 
          representatives (a warning will be raised).
        * ([#191](https://github.com/Ecogenomics/GTDBTk/issues/191)) Added a new command `ani_rep` which calculates the ANI of input genomes to all GTDB 
          representative genomes. 
            * This command uses [Mash](https://github.com/marbl/Mash) in a pre-filtering step. If pre-filtering is enabled (default)
            then `mash` will need to be on the system path. To disable pre-filtering use the `--no_mash` flag.
        * ([#230](https://github.com/Ecogenomics/GTDBTk/issues/235)) Improved how markers are used in determining the correct domain, and gene selection for the alignment.

**Note (Dec 12, 2019)**:
* GTDB-Tk v1.0.2 has been released
    * Fixed an issue where FastANI threads would timeout with `FastANI returned a non-zero exit code.`
        * Versions affected: 1.0.0, and 1.0.1.
    
**Note (Dec 5, 2019)**:
* GTDB-Tk v1.0.0 has been released
    * Migrated to **Python 3**, you must be running at least **Python 3.6** or later to use this version.
    * `check_install` now does an exhaustive check of the reference data.
    * Resolved an issue where gene calling would fail for low quality genomes (#192).
    * Improved FastANI multiprocessing performance.
    * Third party software versions are reported where possible.

**Note (Nov 15, 2019)**:
* GTDB-Tk v0.3.3 has been released
    * A bug has been fixed which affected `classify` and `classify_wf` when using the `--batchfile`
      argument with genome IDs that differed from the FASTA filename. This issue resulted in 
      the assigned taxonomy being derived only from tree placement without any ANI 
      calculations being considered. Consequently, in some cases genomes may have been classified as a new
      species within a genus when they should have been assigned to an existing species. If you have genomes
      with species assignments this bug did not impact you.
    * Progress is now displayed for: hmmalign, and pplacer.
    * Fixed an issue where the `root` command could not be run independently.
    * Improved MSA masking performance.

**Note (Nov 08, 2019)**:
* Python 2 is reaching [end of life](https://pythonclock.org/) on January 1, 2020
    * GTDB-Tk will be ported to Python 3 to accommodate deprecation of Python 2. 
    * The official GTDB-Tk release will be updated to v1.0.0 on Dec 01, 2019.
    * GTDB-Tk v1.0.0 will require Python 3 and there are no plans to support Python 2 moving forward. Apologies for any issues this may cause.
  
**Note (July 12, 2019)**:
* GTDB-Tk v0.3.2 has been released
    * FastANI calculations are more robust.
    * Optimisation of RED calculations.
    * Improved output messages when errors are encountered.

**Note (July 08, 2019)**:
* GTDB-Tk v0.3.1 has been released:
  * Pplacer taxonomy is now available in the summary file.
  * FastANI species assignment will be selected over phylogenetic placement (Topology case).
  
**Note (June 21, 2019)**:
* GTDB-Tk v0.3.0 has been released:
  * Best translation table displayed in summary file.
  * GTDB-Tk now supports gzipped genomes as inputs (--extension .gz).
  * By default, GTDB-Tk uses precalculated RED values.
  * New option to recalculate RED value during classify step (--recalculate_red).
  * New option to export the untrimmed reference MSA files.
  * New option to skip_trimming during align step.
  * New option to use a custom taxonomy file when rooting a tree.
  * New [FAQ](docs/faq.md) page available.
  * New output structure.
  * This version requires a new version of the GTDB-Tk data package (gtdbtk_r89_data.tar.gz) available [here](https://data.ace.uq.edu.au/public/gtdb/data/releases/release89/89.0/)
  
**Note (April 12, 2019)**:
- A GTDB-Tk reference data package (gtdbtk.r86_v3_data.tar.gz) is available [here](https://data.ace.uq.edu.au/public/gtdbtk/release_86/).
- This package contains an additional set of genomes based on GTDB R86.

**Note (March 03, 2019)**:
* GTDB-Tk v0.2.1 has been released:
  * Species classification is now based strictly on the ANI to reference genomes
  * The "classify" function now reports the closest reference genome in the summary file even if the ANI is <95%
  * The summary.tsv file has 4 new columns: aa_percent, red_values, fastani_reference_radius, and warnings
  * By default, the "align" function now performs the same MSA trimming used by the GTDB
  * New pplacer support for writing to a scratch file (--mmap-file option)
  * Random seed option for MSA trimming has been added to allow for reproducible results 
  * Configuration of the data directory is now set using the environmental variable GTDBTK_DATA_PATH (see [pip installation](#pip-installation))
  * Perl dependencies has been removed
  * Python libraries biolib, mpld3 and jinja have been removed
  * This version requires a new version of the GTDB-Tk data package (gtdbtk.r86_v2_data.tar.gz) available [here](https://data.ace.uq.edu.au/public/gtdbtk/release_86/)  

**Note (Sept. 20, 2018)**:
- GTDB-Tk v0.1.3 has been released and addresses an issue with species assignments based on the placement of genomes in the reference tree. This impacted species assignment when submitting multiple closely related genomes. Species assignments reported by ANI were not impacted.

**Note (Aug. 30, 2018)**:
- A new version of the data (release 86) is available under [this link](https://data.ace.uq.edu.au/public/gtdbtk/release_86/).
- This new version is required to run GTDB-Tk v0.1.0+
