.. _installing/bioconda:

Bioconda
========

Step 1: Install anaconda (if not already done)
----------------------------------------------

Ensure that ``conda`` on the system path. It is recommended to download `miniconda <https://docs.conda.io/en/latest/miniconda.html>`_.


Step 2: Create the GTDB-Tk environment
--------------------------------------

.. note:: It is strongly recommended to create a new GTDB-Tk environment for each version of GTDB-Tk released.

GTDB-Tk package requires third-party packages from the ``conda-forge`` and ``bioconda`` channels.


.. code-block:: bash

    # latest version
    conda create -n gtdbtk -c conda-forge -c bioconda gtdbtk

    # specific version (replace 1.3.0 with the version you wish to install, recommended)
    conda create -n gtdbtk-1.3.0 -c conda-forge -c bioconda gtdbtk=1.3.0


Step 3: Download and alias the GTDB-Tk reference data
-----------------------------------------------------

GTDB-Tk requires an environment variable named ``GTDBTK_DATA_PATH`` to be set to the directory
containing the unarchived :ref:`installing#gtdbtk-reference-data`.

Automatically
^^^^^^^^^^^^^

The conda package is bundled with a script ``download-db.sh`` `(located here) <https://github.com/bioconda/bioconda-recipes/blob/master/recipes/gtdbtk/download-db.sh>`_
that will automatically download, and extract the GTDB-Tk reference data. The script will be on the system path so simply run:

.. code-block:: bash

    download-db.sh



Manually
^^^^^^^^

You can automatically alias ``GTDBTK_DATA_PATH`` whenever the environment is activated by editing ``{gtdbtk environment path}/etc/conda/activate.d/gtdbtk.sh``, e.g.:

.. code-block:: bash

    # Determine the GTDB-Tk environment path
    conda activate gtdbtk-1.3.0
    which gtdbtk
    # /miniconda3/envs/gtdbtk-1.3.0/bin/gtdbtk

    # Edit the activate file
    echo "export GTDBTK_DATA_PATH=/path/to/release/package/" > /miniconda3/envs/gtdbtk-1.3.0/etc/conda/activate.d/gtdbtk.sh
