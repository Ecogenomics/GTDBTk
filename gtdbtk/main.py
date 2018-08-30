###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

import os
import logging
import sys

from markers import Markers
from classify import Classify
from misc import Misc
from reroot_tree import RerootTree
import config.config as Config

from biolib.common import (check_dir_exists,
                           check_file_exists,
                           make_sure_path_exists,
                           remove_extension)
from biolib.taxonomy import Taxonomy
from biolib.external.execute import check_dependencies
from biolib.external.fasttree import FastTree


class OptionsParser():

    def __init__(self, version):
        """Initialization."""

        self.version = version

        self.logger = logging.getLogger('timestamp')

        self.logger.warning(
            "Results are still being validated and taxonomic assignments may be incorrect! Use at your own risk!")

    def _verify_genome_id(self, genome_id):
        """Ensure genome ID will be valid in Newick tree."""

        invalid_chars = set('()[],;=')
        if any((c in invalid_chars) for c in genome_id):
            self.logger.error('Invalid genome ID: %s' % genome_id)
            self.logger.error(
                'The following characters are invalid: %s' % ' '.join(invalid_chars))
            sys.exit(-1)

        return True

    def _genomes_to_process(self, genome_dir, batchfile, extension):
        """Get genomes to process.

        Parameters
        ----------
        genome_dir : str
          Directory containing genomes.
        batchfile : str
          File describing genomes.
        extension : str
          Extension of files to process.

        Returns
        -------
        genomic_files : d[genome_id] -> FASTA file
            Map of genomes to their genomic FASTA files.
        """

        genomic_files = {}
        if genome_dir:
            for f in os.listdir(genome_dir):
                if f.endswith(extension):
                    genome_id = remove_extension(f)
                    genomic_files[genome_id] = os.path.join(genome_dir, f)

        elif batchfile:
            for line_no, line in enumerate(open(batchfile, "rb")):
                line_split = line.strip().split("\t")
                if line_split[0] == '':
                    continue  # blank line

                if len(line_split) != 2:
                    self.logger.error(
                        'Batch file must contain exactly 2 columns.')
                    sys.exit(-1)

                genome_file, genome_id = line_split
                self._verify_genome_id(genome_id)

                if genome_file is None or genome_file == '':
                    self.logger.error(
                        'Missing genome file on line %d.' % line_no + 1)
                    self.exit(-1)
                elif genome_id is None or genome_id == '':
                    self.logger.error(
                        'Missing genome ID on line %d.' % line_no + 1)
                    self.exit(-1)
                elif genome_id in genomic_files:
                    self.logger.error(
                        'Genome ID %s appear multiple times.' % genome_id)
                    self.exit(-1)

                genomic_files[genome_id] = genome_file

        for genome_key in genomic_files.iterkeys():
            if genome_key.startswith("RS_") or genome_key.startswith("GB_") or genome_key.startswith("UBA"):
                self.logger.error(
                    "Submitted genomes start with the same prefix (RS_,GB_,UBA) as reference genomes in GTDB-Tk. This will cause issues for downstream analysis.")
                sys.exit(-1)

        if len(genomic_files) == 0:
            if genome_dir:
                self.logger.warning(
                    'No genomes found in directory: %s. Check the --extension flag used to identify genomes.' % genome_dir)
            else:
                self.logger.warning(
                    'No genomes found in batch file: %s. Please check the format of this file.' % batchfile)
            sys.exit(-1)

        return genomic_files

    def _marker_set_id(self, bac120_ms, ar122_ms, rps23_ms):
        """Get unique identifier for marker set."""

        if bac120_ms:
            marker_set_id = "bac120"
        elif ar122_ms:
            marker_set_id = "ar122"
        elif rps23_ms:
            marker_set_id = "rps23"

        return marker_set_id

    def identify(self, options):
        """Identify marker genes in genomes."""

        if options.genome_dir:
            check_dir_exists(options.genome_dir)

        if options.batchfile:
            check_file_exists(options.batchfile)

        make_sure_path_exists(options.out_dir)

        genomes = self._genomes_to_process(
            options.genome_dir, options.batchfile, options.extension)

        markers = Markers(options.cpus)
        markers.identify(genomes,
                         options.out_dir,
                         options.prefix)

        self.logger.info('Done.')

    def align(self, options):
        """Create MSA from marker genes."""

        check_dir_exists(options.identify_dir)
        make_sure_path_exists(options.out_dir)

        if not hasattr(options, 'outgroup_taxon'):
            options.outgroup_taxon = None

        markers = Markers(options.cpus)
        markers.align(options.identify_dir,
                      options.taxa_filter,
                      options.min_perc_aa,
                      options.custom_msa_filters,
                      options.consensus,
                      options.min_perc_taxa,
                      options.out_dir,
                      options.prefix,
                      options.outgroup_taxon)

        self.logger.info('Done.')

    def infer(self, options):
        """Infer tree from MSA."""

        self.logger.warning("Tree inference is still under development!")

        check_file_exists(options.msa_file)
        make_sure_path_exists(options.out_dir)

        if (options.cpus > 1):
            check_dependencies(['FastTreeMP'])
        else:
            check_dependencies(['FastTree'])

        self.logger.info(
            'Inferring tree with FastTree using %s+GAMMA.' % options.prot_model)
        fasttree = FastTree(multithreaded=(options.cpus > 1))

        tree_unrooted_output = os.path.join(
            options.out_dir, options.prefix + options.suffix + '.unrooted.tree')
        tree_log = os.path.join(options.out_dir, options.prefix + '.tree.log')
        tree_output_log = os.path.join(options.out_dir, 'fasttree.log')
        fasttree.run(options.msa_file,
                     'prot',
                     options.prot_model,
                     tree_unrooted_output,
                     tree_log,
                     tree_output_log)

        self.logger.info('Done.')

    def classify(self, options):
        """Determine taxonomic classification of genomes."""

        check_dir_exists(options.align_dir)
        make_sure_path_exists(options.out_dir)

        genomes = self._genomes_to_process(
            options.genome_dir, options.batchfile, options.extension)

        classify = Classify(options.cpus)
        classify.run(genomes,
                     options.align_dir,
                     options.out_dir,
                     options.prefix,
                     options.debug)

        self.logger.info('Done.')

    def trim_msa(self, options):
        """ Trim an untrimmed archaea or bacterial MSA file."""
        if options.reference_mask in ['bac', 'arc']:
            mask_type = "reference"
            mask_id = options.reference_mask
        else:
            mask_type = "file"
            mask_id = options.mask_file
        misc = Misc()
        misc.trim_msa(options.untrimmed_msa, mask_type,
                      mask_id, options.output)
        self.logger.info('Done.')

    def root(self, options):
        """Root tree using outgroup."""
        self.logger.warning("Tree rooting is still under development!")

        check_file_exists(options.input_tree)

        gtdb_taxonomy = Taxonomy().read(Config.TAXONOMY_FILE)

        self.logger.info('Identifying genomes from the specified outgroup.')
        outgroup = set()
        for genome_id, taxa in gtdb_taxonomy.iteritems():
            if options.outgroup_taxon in taxa:
                outgroup.add(genome_id)

        reroot = RerootTree()
        reroot.root_with_outgroup(options.input_tree,
                                  options.output_tree,
                                  outgroup)

        self.logger.info('Done.')

    def check_install(self):
        """ Verify is all gtdb data files are present to run GTDB-Tk."""
        self.logger.warning("Running install verification")
        misc = Misc()
        misc.check_install()
        self.logger.info('Done.')

    def decorate(self, options):
        """Decorate tree with GTDB taxonomy."""

        check_file_exists(options.input_tree)

        # Config.TAXONOMY_FILE
        self.logger.warning('NOT YET IMPLEMENTED!')
        self.logger.info('Done.')

    def parse_options(self, options):
        """Parse user options and call the correct pipeline(s)"""

        if(options.subparser_name == 'de_novo_wf'):
            check_dependencies(['prodigal', 'hmmalign'])
            if (options.cpus > 1):
                check_dependencies(['FastTreeMP'])
            else:
                check_dependencies(['FastTree'])

            self.identify(options)

            options.identify_dir = options.out_dir
            self.align(options)

            if options.bac120_ms:
                options.suffix = ".bac120"
            else:
                options.suffix = ".ar122"

            options.msa_file = os.path.join(options.out_dir,
                                            options.prefix + options.suffix + ".msa.fasta")
            self.infer(options)

            options.input_tree = os.path.join(options.out_dir,
                                              options.prefix + options.suffix + ".unrooted.tree")
            options.output_tree = os.path.join(options.out_dir,
                                               options.prefix + options.suffix + ".rooted.tree")
            self.root(options)

            self.decorate(options)
        elif(options.subparser_name == 'classify_wf'):
            check_dependencies(
                ['prodigal', 'hmmalign', 'pplacer', 'guppy', 'fastANI'])
            self.identify(options)

            options.identify_dir = options.out_dir
            options.align_dir = options.out_dir
            options.taxa_filter = None
            options.custom_msa_filters = False
            options.consensus = None
            options.min_perc_taxa = None
            self.align(options)

            self.classify(options)
        elif (options.subparser_name == 'identify'):
            self.identify(options)
        elif(options.subparser_name == 'align'):
            self.align(options)
        elif(options.subparser_name == 'infer'):
            self.infer(options)
        elif(options.subparser_name == 'classify'):
            self.classify(options)
        elif(options.subparser_name == 'root'):
            self.root(options)
        elif(options.subparser_name == 'decorate'):
            self.decorate(options)
        elif(options.subparser_name == 'trim_msa'):
            self.trim_msa(options)
        elif(options.subparser_name == 'check_install'):
            self.check_install()
        else:
            self.logger.error('Unknown GTDB-Tk command: "' +
                              options.subparser_name + '"\n')
            sys.exit()

        return 0
