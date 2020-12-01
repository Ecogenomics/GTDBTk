.. _installing/pip:

pip
===


Step 1: Install third-party dependencies
----------------------------------------

Ensure that the software described in :ref:`installing#third-party-software` are on the system path.


Step 2: Install GTDB-Tk via pip
-------------------------------

Note: It is strongly recommended to create a `virtual environment <https://packaging.python.org/guides/installing-using-pip-and-virtual-environments/>`_
for each version of GTDB-Tk installed.


Once the third-party dependencies have been installed, install GTDB-Tk via `pip <https://pypi.python.org/pypi/gtdbtk>`_:

.. code-block:: bash

    python -m pip install gtdbtk


Step 3: Download and alias the GTDB-Tk reference data
-----------------------------------------------------

GTDB-Tk requires an environment variable named ``GTDBTK_DATA_PATH`` to be set to the directory
containing the unarchived :ref:`installing#gtdbtk-reference-data`.

.. code-block:: bash

    export GTDBTK_DATA_PATH=/path/to/release/package/

You can permanently save this variable as described `here <https://unix.stackexchange.com/questions/26047/how-to-correctly-add-a-path-to-path>`_.
