*******
GTDB-Tk
*******


GTDB-Tk is a software toolkit for assigning objective taxonomic classifications to bacterial and
archaeal genomes based on the Genome Database Taxonomy `GTDB <https://gtdb.ecogenomic.org/>`_.
It is designed to work with recent advances that allow hundreds or thousands of metagenome-assembled
genomes (MAGs) to be obtained directly from environmental samples. It can also be applied to isolate
and single-cell genomes. The GTDB-Tk is open source and released under the
`GNU General Public License (Version 3) <https://www.gnu.org/licenses/gpl-3.0.en.html>`_.

Notifications about GTDB-Tk releases will be available through the GTDB Twitter account
`<https://twitter.com/ace_gtdb>`_.

Please post questions and issues related to GTDB-Tk on the Issues section of the GitHub repository.
Questions related to the `GTDB <https://gtdb.ecogenomic.org/>`_ should be sent to the
`GTDB team <https://gtdb.ecogenomic.org/about>`_.


Running GTDB-Tk
===============

1. Install GTDB-Tk (or use the third-party web application) (:ref:`installing`)

2. Access the help documentation :ref:`commands`, or view the program help menu: ``gtdbtk -h``

* Note: Individual help can be accessed via the specific command, e.g.: ``gtdbtk classify_wf -h``


Citing GTDB-Tk
==============
We encourage you to cite GTDB-Tk and the third-party dependencies as described in :ref:`references`.


.. toctree::
   :caption: Getting started
   :maxdepth: 1

   announcements
   installing/index
   faq


.. toctree::
   :caption: Running GTDB-Tk
   :maxdepth: 1

   commands/index
   files/index
   examples/classify_wf


.. toctree::
   :caption: About
   :maxdepth: 1

   changelog
   references
