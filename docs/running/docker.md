---
layout: default
title: Docker
nav_order: 4
parent: Running GTDB-Tk
has_children: false
---

# Docker

[![Docker Image Version (tag latest semver)](https://img.shields.io/docker/v/ecogenomic/gtdbtk/latest?color=299bec&label=docker)](https://hub.docker.com/repository/docker/ecogenomic/gtdbtk)
[![Docker Pulls](https://img.shields.io/docker/pulls/ecogenomic/gtdbtk?color=299bec&label=pulls)](https://hub.docker.com/repository/docker/ecogenomic/gtdbtk)

## Prerequisites

The docker image requires two directories to be mounted.

### 1. GTDB-Tk package data:
Within the container, GTDB-Tk will look for the reference data under 

```bash
/refdata
```

Therefore, you will need to specify where you have downloaded the reference data on your local system, e.g.:

```bash
-v /host/release89:/refdata
```

### 2. Input/output data:
In order for the container to access any input genomes, or for you to be able to access any results, 
you will need to mount a local directory for the container to use. Within the container, GTDB-Tk will be able to write to:
 
 ```bash
/data
 ```

Therefore, you will need to specify where you are willing to allow GTDB-Tk to read/write on your local system, e.g.:

```bash
-v /host/gtdbtk_output:/data
 ```

## Running the container

Putting this all together, a command may look like (note that since `gtdbtk` is the entrypoint, this can be ommited):

```bash
docker run -v /host/gtdbtk_output:/data -v /host/release89:/refdata ecogenomic/gtdbtk --help
```

or

```bash
docker run -v /host/gtdbtk_output:/data -v /host/release89:/refdata ecogenomic/gtdbtk test --out_dir /data/output
```
