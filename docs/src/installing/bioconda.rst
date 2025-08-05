.. _installing/bioconda:

Bioconda
========

Step 1: Install conda (if not already done)
-------------------------------------------

We strongly recommend using `Mamba <https://mamba.readthedocs.io/en/latest/installation.html>`_ (much faster!) over `miniconda <https://docs.conda.io/en/latest/miniconda.html>`_/`conda <https://docs.conda.io/projects/conda/en/stable/user-guide/install/index.html>`_, but all will work.


Step 2: Create the GTDB-Tk environment
--------------------------------------

It is strongly recommended to create a new conda environment for each version of GTDB-Tk released.

.. warning:: You must always specify the version of GTDB-Tk, as conda may try to install a **very old version (v1.0.2)**.


GTDB-Tk requires third-party packages from the ``conda-forge`` and ``bioconda`` channels, make sure to
specify those channels in that order!

.. code-block:: bash

    # NOTE: replace |release| with the version you wish to install

    # using conda
    conda create -n gtdbtk-|release| -c conda-forge -c bioconda gtdbtk=|release|

    # using mamba (alternative)
    mamba create -n gtdbtk-|release| -c conda-forge -c bioconda gtdbtk=|release|

Step 3: Download and alias the GTDB-Tk reference data
-----------------------------------------------------

GTDB-Tk requires an environment variable named ``GTDBTK_DATA_PATH`` to be set to the directory
containing the unarchived :ref:`installing#gtdbtk-reference-data`.

Automatically
^^^^^^^^^^^^^

The conda package is bundled with a script ``download-db.sh`` `(source) <https://github.com/bioconda/bioconda-recipes/blob/master/recipes/gtdbtk/download-db.sh>`_
that will automatically download, and extract the GTDB-Tk reference data. The script will be on the system path so simply run:

.. code-block:: bash

    download-db.sh



Manually
^^^^^^^^

You can automatically alias ``GTDBTK_DATA_PATH`` whenever the environment is activated by
`setting environment-specific variables <https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#setting-environment-variables>`_, e.g.:

.. code-block:: bash

    # Activate the GTDB-Tk conda environment
    conda activate gtdbtk-|release|

    # Set the environment variable to the directory containing the GTDB-Tk reference data
    conda env config vars set GTDBTK_DATA_PATH="/path/to/unarchived/gtdbtk/data";
