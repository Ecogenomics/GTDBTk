---
layout: default
title: pip
nav_order: 2
parent: Running GTDB-Tk
has_children: false
---

# pip

Python >=3.6
{: .label }

[![PyPI](https://img.shields.io/pypi/v/gtdbtk.svg)](https://pypi.python.org/pypi/gtdbtk)
[![PyPI Downloads](https://pepy.tech/badge/gtdbtk)](https://pepy.tech/project/gtdbtk)

## Installation

Once the third-party dependencies have been installed, GTDB-Tk can be installed using [pip](https://pypi.python.org/pypi/gtdbtk):
```bash
python -m pip install gtdbtk
```

## Post-install

GTDB-Tk requires an environment variable named `GTDBTK_DATA_PATH` to be set to the directory 
containing the unarchived GTDB-Tk reference data.
```
export GTDBTK_DATA_PATH=/path/to/release/package/
```
You can permanently save this variable as described [here](https://unix.stackexchange.com/questions/26047/how-to-correctly-add-a-path-to-path).
