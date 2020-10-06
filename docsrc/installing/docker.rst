.. _installing/docker:

Docker
======


Step 1: Install Docker (if not already done)
----------------------------------------------

Ensure that ``docker`` on the system path, see `Docker <https://www.docker.com/get-started>`_ for more information.


Step 2: Download the GTDB-Tk reference data
-------------------------------------------

Download an unarchive the :ref:`installing#gtdbtk-reference-data` to a directory of your choice, e.g: ``/host/release_data/``


Step 3: Create a local I/O directory for GTDB-Tk
------------------------------------------------

You will need to create a directory on the host machine that the container can access. The container
will use this to read input genomes, and write output files, e.g.: ``/host/gtdbtk_io/``

You should either copy or symlink your genomes to this directory, e.g. you may have the following contents:

.. code-block:: text

     /host/gtdbtk_io/genomes/genome_a.fna
     /host/gtdbtk_io/genomes/genome_b.fna


Step 4: Run the container
-------------------------

The two host directories created in Step 2, and Step 3 will be mapped to the container via the ``-v`` flag:

* ``/host/release_data:/refdata``
* ``/host/gtdbtk_io:/data``

For example, the classify workflow for genome_a, and genome_b an be run as follows:

.. code-block:: bash

     docker run -v /host/gtdbtk_io:/data -v /host/release_data:/refdata ecogenomic/gtdbtk classify_wf --genome_dir /data/genomes --out_dir /data/output

Note:

* When referring to ``/host/gtdbtk_io`` in the GTDB-Tk arguments, it has been aliased to ``/data`` as specified by the ``-v`` flag.
* The alias ``/data`` and ``/refdata`` cannot be changed.

