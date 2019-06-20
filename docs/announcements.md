
## Previous announcements


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