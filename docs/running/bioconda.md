---
layout: default
title: Bioconda
nav_order: 3
parent: Running GTDB-Tk
has_children: false
---

# Bioconda

[![Bioconda](https://img.shields.io/conda/vn/bioconda/gtdbtk.svg?color=43b02a)](https://anaconda.org/bioconda/gtdbtk)
[![BioConda Downloads](https://img.shields.io/conda/dn/bioconda/gtdbtk.svg?style=flag&label=downloads&color=43b02a)](https://anaconda.org/bioconda/gtdbtk)

## Installation

1. Create the conda environment:
    * Latest version: `conda create -n gtdbtk -c bioconda gtdbtk`
    * Specific version: `conda create -n gtdbtk-1.1.0 -c bioconda gtdbtk=1.1.0`
2. Activate the environment: `conda activate gtdbtk`

## Post-install

1. Download the reference package either [manually](https://github.com/Ecogenomics/GTDBTk#gtdb-tk-reference-data) 
or by running `download-db.sh`.
2. Set the `GTDBTK_DATA_PATH` environment variable in 
`{gtdbtk environment path}/etc/conda/activate.d/gtdbtk.sh` to the reference package location.

