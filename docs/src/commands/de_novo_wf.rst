.. _commands/de_novo_wf:

de_novo_wf
==========

For arguments and output files, see each of the individual steps:

* :ref:`commands/identify`
* :ref:`commands/align`
* :ref:`commands/infer`
* :ref:`commands/root`
* :ref:`commands/decorate`


The *de novo* workflow infers new bacterial and archaeal trees containing all user supplied and GTDB-Tk reference genomes.
The classify workflow is recommended for obtaining taxonomic classifications, and this workflow only recommended if
a *de novo* domain-specific trees are desired. **One should take the taxonomic assignments as a guide, but not as final classifications**. In particular, no effort is made to resolve the taxonomic assignment of lineages composed exclusively of user submitted genomes.

This workflow consists of five steps: ``identify``, ``align``, ``infer``, ``root``,
and ``decorate``.

The ``identify`` and ``align`` steps are the same as in the classify workflow.

The ``infer`` step uses `FastTree <http://www.microbesonline.org/fasttree/>`_ with the WAG+GAMMA models to calculate independent, *de novo* bacterial and archaeal trees.
These trees can then be rooted using a user specified outgroup and decorated with the GTDB taxonomy.

The *de novo* workflow can be run as follows:

.. code-block:: bash

    gtdbtk de_novo_wf --genome_dir <my_genomes> --<marker_set> --outgroup_taxon <outgroup> --out_dir <output_dir>


This will process all genomes in <my_genomes> using the specified marker set and place the results in <output_dir>.
Only genomes previously identified as being bacterial (archaeal) should be included when using the bacterial (archaeal) marker set.
The tree will be rooted with the <outgroup> taxon (typically a phylum in the domain-specific tree) as required for
correct decoration of the tree. In general, we suggest the resulting tree be treated as unrooted when interpreting results.
Identical to the classify workflow, the location of genomes can also be specified using a batch file with the ``--batchfile`` flag.


The workflow supports several optional flags, including:

* cpus: maximum number of CPUs to use
* min_perc_aa: filter genomes with an insufficient percentage of AA in the MSA (default: 50)
* taxa_filter: filter genomes to taxa within specific taxonomic groups
* prot_model: protein substitution model for tree inference (LG or WAG; default: WAG)

For other flags please consult the command line interface.


Arguments
---------

.. argparse::
   :module: gtdbtk.cli
   :func: get_main_parser
   :prog: gtdbtk
   :path: de_novo_wf
   :nodefaultconst:



Example
-------

Input
^^^^^

.. code-block:: bash

    gtdbtk de_novo_wf --genome_dir genomes/ --outgroup_taxon p__Undinarchaeota --archaea --out_dir de_novo_wf --cpus 3
    
    gtdbtk de_novo_wf --genome_dir genomes/ --outgroup_taxon p__Chloroflexota --bacteria  --taxa_filter p__Firmicutes --out_dir de_novo_output

    #Skip GTDB reference genomes ( requires --custom_taxonomy_file for outgrouping)
    gtdbtk de_novo_wf --genome_dir genomes/ --outgroup_taxon p__Customphylum --bacteria --custom_taxonomy_file custom_taxonomy.tsv --out_dir de_novo_output

    #Use a subset of GTDB reference genomes (p__Firmicutes) and outgroup on a custom Phylum (p__Customphylum)
    gtdbtk de_novo_wf --genome_dir genomes/ --taxa_filter p__Firmicutes --outgroup_taxon p__Customphylum --bacteria --custom_taxonomy_file custom_taxonomy.tsv --out_dir de_novo_output

Custom Taxonomy Format
^^^^^^^^^^^^^^^^^^^^^^
The custom taxonomy file is a Tab-delimited file with the first column listing user genomes (i.e Fasta filename without the extension)
and the second column listing the standardized 7-rank taxonomy.

.. code-block:: bash

    #For genome_1.fna, genome_2.fna and genome_3.fna
    genome_1    d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Salmonella;s__Salmonella enterica
    genome_2    d__Bacteria;p__Actinobacteriota;c__Actinomycetia;o__Mycobacteriales;f__Mycobacteriaceae;g__Mycobacterium;s__Mycobacterium tuberculosis
    genome_3    d__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Streptococcaceae;g__Streptococcus;s__Streptococcus pyogenes
