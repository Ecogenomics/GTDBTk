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

from __future__ import print_function

import os
import logging
import sys
import shutil

from gtdbtk.tools import symlink_f
from markers import Markers
from classify import Classify
from misc import Misc
from reroot_tree import RerootTree
import config.config as Config
from gtdbtk.config.output import *

from biolib_lite.common import (check_dir_exists,
                                check_file_exists,
                                make_sure_path_exists,
                                remove_extension)
from biolib_lite.taxonomy import Taxonomy
from biolib_lite.execute import check_dependencies
from gtdbtk.exceptions import *


class OptionsParser(object):

    def __init__(self, version):
        """Initialization."""
        self.logger = logging.getLogger('timestamp')
        self.version = version

        self.logger.warning(
            "Results are still being validated and taxonomic assignments may be incorrect! Use at your own risk!")

        self._check_package_compatibility()

    def _check_package_compatibility(self):
        """Check that GTDB-Tk is using the most up-to-date reference package."""
        pkg_ver = float(Config.VERSION_DATA.replace('r', ''))
        min_ver = float(Config.MIN_REF_DATA_VERSION.replace('r', ''))
        self.logger.info('Using GTDB-Tk reference data version {}: {}'
                         .format(Config.VERSION_DATA, Config.GENERIC_PATH))
        if pkg_ver < min_ver:
            self.logger.warning('You are not using the reference data intended '
                                'for this release: {}'.format(Config.MIN_REF_DATA_VERSION))

    def _verify_genome_id(self, genome_id):
        """Ensure genome ID will be valid in Newick tree."""

        invalid_chars = set('()[],;=')
        if any((c in invalid_chars) for c in genome_id):
            self.logger.error('Invalid genome ID: %s' % genome_id)
            self.logger.error('The following characters are invalid: %s' % ' '.join(invalid_chars))
            raise GenomeNameInvalid
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
                    self.logger.error('Batch file must contain exactly 2 columns.')
                    raise GenomeBatchfileMalformed

                genome_file, genome_id = line_split
                self._verify_genome_id(genome_id)

                if genome_file is None or genome_file == '':
                    self.logger.error('Missing genome file on line %d.' % (line_no + 1))
                    raise GenomeBatchfileMalformed
                elif genome_id is None or genome_id == '':
                    self.logger.error('Missing genome ID on line %d.' % (line_no + 1))
                    raise GenomeBatchfileMalformed
                elif genome_id in genomic_files:
                    self.logger.error('Genome ID %s appear multiple times.' % genome_id)
                    raise GenomeBatchfileMalformed
                if genome_file in genomic_files.values():
                    self.logger.warning('Genome file appears multiple times: %s' % genome_file)

                genomic_files[genome_id] = genome_file

        for genome_key in genomic_files.iterkeys():
            if genome_key.startswith("RS_") or genome_key.startswith("GB_") or genome_key.startswith("UBA"):
                self.logger.error(
                    "Submitted genomes start with the same prefix (RS_,GB_,UBA) as reference genomes in GTDB-Tk. This will cause issues for downstream analysis.")
                raise GenomeNameInvalid

        if len(genomic_files) == 0:
            if genome_dir:
                self.logger.error(
                    'No genomes found in directory: %s. Check the --extension flag used to identify genomes.' % genome_dir)
            else:
                self.logger.error(
                    'No genomes found in batch file: %s. Please check the format of this file.' % batchfile)
            raise NoGenomesFound

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
                         options.prefix,
                         options.force)

        self.logger.info('Done.')

    def align(self, options):
        """Create MSA from marker genes."""

        check_dir_exists(options.identify_dir)
        make_sure_path_exists(options.out_dir)

        if not hasattr(options, 'outgroup_taxon'):
            options.outgroup_taxon = None

        markers = Markers(options.cpus)
        markers.align(options.identify_dir,
                      options.skip_gtdb_refs,
                      options.taxa_filter,
                      options.min_perc_aa,
                      options.custom_msa_filters,
                      options.skip_trimming,
                      options.rnd_seed,
                      options.cols_per_gene,
                      options.min_consensus,
                      options.max_consensus,
                      options.min_perc_taxa,
                      options.out_dir,
                      options.prefix,
                      options.outgroup_taxon)

        self.logger.info('Done.')

    def infer(self, options):
        """Infer tree from MSA."""

        check_file_exists(options.msa_file)
        make_sure_path_exists(options.out_dir)

        if options.cpus > 1:
            check_dependencies(['FastTreeMP'])
            os.environ['OMP_NUM_THREADS'] = '%d' % options.cpus
        else:
            check_dependencies(['FastTree'])

        if hasattr(options, 'suffix'):
            output_tree = os.path.join(options.out_dir,
                                       PATH_MARKER_UNROOTED_TREE.format(prefix=options.prefix, marker=options.suffix))
            tree_log = os.path.join(options.out_dir,
                                    PATH_MARKER_TREE_LOG.format(prefix=options.prefix, marker=options.suffix))
            fasttree_log = os.path.join(options.out_dir,
                                        PATH_MARKER_FASTTREE_LOG.format(prefix=options.prefix, marker=options.suffix))
        else:
            output_tree = os.path.join(options.out_dir, PATH_UNROOTED_TREE.format(prefix=options.prefix))
            tree_log = os.path.join(options.out_dir, PATH_TREE_LOG.format(prefix=options.prefix))
            fasttree_log = os.path.join(options.out_dir, PATH_FASTTREE_LOG.format(prefix=options.prefix))

        make_sure_path_exists(os.path.dirname(output_tree))
        make_sure_path_exists(os.path.dirname(tree_log))
        make_sure_path_exists(os.path.dirname(fasttree_log))

        if options.prot_model == 'JTT':
            model_str = ''
        elif options.prot_model == 'WAG':
            model_str = ' -wag'
        elif options.prot_model == 'LG':
            model_str = ' -lg'

        support_str = ''
        if options.no_support:
            support_str = ' -nosupport'

        gamma_str = ' -gamma'
        gamma_str_info = '+GAMMA'
        if options.no_gamma:
            gamma_str = ''
            gamma_str_info = ''

        self.logger.info(
            'Inferring tree with FastTree using {}.'.format(options.prot_model, gamma_str_info))

        cmd = '-quiet%s%s%s -log %s %s > %s 2> %s' % (support_str,
                                                      model_str,
                                                      gamma_str,
                                                      tree_log,
                                                      options.msa_file,
                                                      output_tree,
                                                      fasttree_log)
        if options.cpus > 1:
            cmd = 'FastTreeMP ' + cmd
        else:
            cmd = 'FastTree ' + cmd
        os.system(cmd)

        self.logger.info('Done.')

    def run_test(self, options):
        """Run test of classify workflow."""

        make_sure_path_exists(options.out_dir)

        output_dir = os.path.join(options.out_dir, 'output')
        genome_test_dir = os.path.join(options.out_dir, 'genomes')
        if os.path.exists(genome_test_dir):
            self.logger.error('Test directory {} already exists. Test must be run with a new directory.'.format(
                                genome_test_dir))
            sys.exit(-1)

        current_path = os.path.dirname(os.path.realpath(__file__))
        input_dir = os.path.join(current_path, 'tests', 'data', 'genomes')

        shutil.copytree(input_dir, genome_test_dir)

        cmd = 'gtdbtk classify_wf --genome_dir {} --out_dir {} --cpus {}'.format(
            genome_test_dir, output_dir, options.cpus)
        print("Command:")
        print(cmd)
        os.system(cmd)
        summary_file = os.path.join(output_dir, PATH_AR122_SUMMARY_OUT.format(prefix='gtdbtk'))

        if not os.path.exists(summary_file):
            print("{} is missing.\nTest has failed.".format(summary_file))
            sys.exit(-1)

        self.logger.info('Test has successfully finished.')

    def classify(self, options):
        """Determine taxonomic classification of genomes."""

        check_dir_exists(options.align_dir)
        make_sure_path_exists(options.out_dir)
        if options.scratch_dir:
            make_sure_path_exists(options.scratch_dir)

        genomes = self._genomes_to_process(
            options.genome_dir, options.batchfile, options.extension)

        classify = Classify(options.cpus)
        classify.run(genomes,
                     options.align_dir,
                     options.out_dir,
                     options.prefix,
                     options.scratch_dir,
                     options.recalculate_red,
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

    def export_msa(self, options):
        """Export the untrimmed archaeal or bacterial MSA file."""
        misc = Misc()
        misc.export_msa(options.domain, options.output)

        self.logger.info('Done.')

    def root(self, options):
        """Root tree using outgroup."""
        self.logger.warning("Tree rooting is still under development!")

        check_file_exists(options.input_tree)

        if options.custom_taxonomy_file:
            check_file_exists(options.custom_taxonomy_file)
            taxonomy = Taxonomy().read(options.custom_taxonomy_file)
        else:
            taxonomy = Taxonomy().read(Config.TAXONOMY_FILE)

        self.logger.info('Identifying genomes from the specified outgroup.')
        outgroup = set()
        for genome_id, taxa in taxonomy.iteritems():
            if options.outgroup_taxon in taxa:
                outgroup.add(genome_id)

        reroot = RerootTree()
        reroot.root_with_outgroup(options.input_tree,
                                  options.output_tree,
                                  outgroup)

        # Symlink to the tree summary file
        if options.suffix == 'bac120':
            symlink_f(PATH_BAC120_ROOTED_TREE.format(prefix=options.prefix),
                       os.path.join(options.out_dir, os.path.basename(PATH_AR122_ROOTED_TREE.format(prefix=options.prefix))))
        elif options.suffix == 'ar122':
            symlink_f(PATH_AR122_ROOTED_TREE.format(prefix=options.prefix),
                       os.path.join(options.out_dir, os.path.basename(PATH_AR122_ROOTED_TREE.format(prefix=options.prefix))))
        else:
            self.logger.error('There was an error determining the marker set.')
            raise GenomeMarkerSetUnknown

        self.logger.info('Done.')

    def check_install(self):
        """ Verify all GTDB-Tk data files are present."""
        self.logger.info("Running install verification")
        misc = Misc()
        misc.check_install()
        self.logger.info('Done.')

    def decorate(self, options):
        """Decorate tree with GTDB taxonomy."""

        check_file_exists(options.input_tree)

        # Config.TAXONOMY_FILE
        self.logger.warning('DECORATE NOT YET IMPLEMENTED!')
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
            options.skip_trimming = False
            self.align(options)

            if options.bac120_ms:
                options.suffix = "bac120"
            else:
                options.suffix = "ar122"

            if options.skip_gtdb_refs:
                if options.suffix == 'bac120':
                    options.msa_file = os.path.join(options.out_dir, PATH_BAC120_USER_MSA.format(prefix=options.prefix))
                elif options.suffix == 'ar122':
                    options.msa_file = os.path.join(options.out_dir, PATH_AR122_USER_MSA.format(prefix=options.prefix))
                else:
                    self.logger.error('There was an error determining the marker set.')
                    raise GenomeMarkerSetUnknown
            else:
                if options.suffix == 'bac120':
                    options.msa_file = os.path.join(options.out_dir, PATH_BAC120_MSA.format(prefix=options.prefix))
                elif options.suffix == 'ar122':
                    options.msa_file = os.path.join(options.out_dir, PATH_AR122_MSA.format(prefix=options.prefix))
                else:
                    self.logger.error('There was an error determining the marker set.')
                    raise GenomeMarkerSetUnknown

            self.infer(options)

            if options.suffix == 'bac120':
                options.input_tree = os.path.join(options.out_dir, PATH_BAC120_UNROOTED_TREE.format(prefix=options.prefix))
                options.output_tree = os.path.join(options.out_dir, PATH_BAC120_ROOTED_TREE.format(prefix=options.prefix))
            elif options.suffix == 'ar122':
                options.input_tree = os.path.join(options.out_dir, PATH_AR122_UNROOTED_TREE.format(prefix=options.prefix))
                options.output_tree = os.path.join(options.out_dir, PATH_AR122_ROOTED_TREE.format(prefix=options.prefix))
            else:
                self.logger.error('There was an error determining the marker set.')
                raise GenomeMarkerSetUnknown

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
            options.skip_trimming = False  # Added here due to the other mutex argument being include above.
            options.min_consensus = None
            options.min_perc_taxa = None
            options.skip_gtdb_refs = False
            options.cols_per_gene = None
            options.max_consensus = None
            options.rnd_seed = None
            options.skip_trimming = False
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
        elif(options.subparser_name == 'export_msa'):
            self.export_msa(options)
        elif(options.subparser_name == 'test'):
            self.run_test(options)
        elif(options.subparser_name == 'check_install'):
            self.check_install()
        else:
            self.logger.error('Unknown GTDB-Tk command: "' +
                              options.subparser_name + '"\n')
            sys.exit()

        return 0
