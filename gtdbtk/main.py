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

import logging
import os
import shutil
import subprocess
import sys
import tempfile
from typing import Dict, Tuple

from tqdm import tqdm

import gtdbtk.config.config as Config
from gtdbtk.ani_rep import ANIRep
from gtdbtk.biolib_lite.common import (check_dir_exists,
                                       check_file_exists,
                                       make_sure_path_exists,
                                       remove_extension)
from gtdbtk.biolib_lite.execute import check_dependencies
from gtdbtk.biolib_lite.logger import colour
from gtdbtk.biolib_lite.taxonomy import Taxonomy
from gtdbtk.classify import Classify
from gtdbtk.config.output import *
from gtdbtk.decorate import Decorate
from gtdbtk.exceptions import *
from gtdbtk.external.fasttree import FastTree
from gtdbtk.infer_ranks import InferRanks
from gtdbtk.io.batchfile import Batchfile
from gtdbtk.io.classify_summary import ClassifySummaryFileAR122
from gtdbtk.markers import Markers
from gtdbtk.misc import Misc
from gtdbtk.reroot_tree import RerootTree
from gtdbtk.tools import symlink_f, get_reference_ids


class OptionsParser(object):

    def __init__(self, version):
        """Initialization.

        Parameters
        ----------
        version : str
            The current version number (e.g. 0.2.2).
        """
        self.logger = logging.getLogger('timestamp')
        self.warnings = logging.getLogger('warnings')
        self.version = version
        self._check_package_compatibility()

        self.genomes_to_process = None

    def _check_package_compatibility(self):
        """Check that GTDB-Tk is using the most up-to-date reference package."""
        pkg_ver = float(Config.VERSION_DATA.replace('r', ''))
        min_ver = float(Config.MIN_REF_DATA_VERSION.replace('r', ''))
        self.logger.info(f'Using GTDB-Tk reference data version '
                         f'{Config.VERSION_DATA}: {Config.GENERIC_PATH}')
        if pkg_ver < min_ver:
            self.logger.warning(colour(f'You are not using the reference data '
                                       f'intended for this release: {Config.MIN_REF_DATA_VERSION}',
                                       ['bright'], fg='yellow'))

    def _verify_genome_id(self, genome_id):
        """Ensure genome ID will be valid in Newick tree.

        Parameters
        ----------
        genome_id : str
            The string representing the genome identifier.

        Returns
        -------
        bool
            True if the genome identifier is legal.

        Raises
        ------
        GenomeNameInvalid
            If the genome identifier contains illegal characters.
        """

        invalid_chars = set('()[],;=')
        if any((c in invalid_chars) for c in genome_id):
            self.logger.error(f'Invalid genome ID: {genome_id}')
            self.logger.error(f'The following characters are invalid: '
                              f'{" ".join(invalid_chars)}')
            raise GenomeNameInvalid(f'Invalid genome ID: {genome_id}')
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

        genomic_files, tln_tables = dict(), dict()
        if genome_dir:
            for f in os.listdir(genome_dir):
                if f.endswith(extension):
                    genome_id = remove_extension(f, extension)
                    genomic_files[genome_id] = os.path.join(genome_dir, f)

        elif batchfile:
            batchfile_fh = Batchfile(batchfile)
            genomic_files, tln_tables = batchfile_fh.genome_path, batchfile_fh.genome_tln

        # Check that all of the genome IDs are valid.
        for genome_key in genomic_files:
            self._verify_genome_id(genome_key)

        # Check that the prefix is valid and the path exists
        invalid_paths = list()
        for genome_key, genome_path in genomic_files.items():

            if not os.path.isfile(genome_path):
                invalid_paths.append((genome_key, genome_path))

        # Report on any invalid paths
        if len(invalid_paths) > 0:
            self.warnings.info(f'Reading from batchfile: {batchfile}')
            self.warnings.error(f'The following {len(invalid_paths)} genomes '
                                f'have invalid paths specified in the batchfile:')
            for g_path, g_gid in invalid_paths:
                self.warnings.info(f'{g_gid}\t{g_path}')
            raise GTDBTkExit(f'There are {len(invalid_paths)} paths in the '
                             f'batchfile which do not exist, see gtdb.warnings.log')

        if len(genomic_files) == 0:
            if genome_dir:
                self.logger.error('No genomes found in directory: %s. Check '
                                  'the --extension flag used to identify '
                                  'genomes.' % genome_dir)
            else:
                self.logger.error('No genomes found in batch file: %s. Please '
                                  'check the format of this file.' % batchfile)
            raise GTDBTkExit

        invalid_genomes = set(genomic_files.keys()) & set(get_reference_ids())
        if len(invalid_genomes) > 0:
            self.warnings.info(f'The following {len(invalid_genomes)} have the '
                               f'same ID as GTDB-Tk reference genomes:')
            for invalid_genome in sorted(invalid_genomes):
                self.warnings.info(invalid_genome)
            raise GTDBTkExit(f'You have {len(invalid_genomes)} genomes with the '
                             f'same id as GTDB-Tk reference genomes, please '
                             f'rename them. See gtdb.warnings.log.')

        return genomic_files, tln_tables

    def _read_taxonomy_files(self, options) -> Dict[str, Tuple[str, str, str, str, str, str, str]]:
        """Read and merge taxonomy files."""

        self.logger.info('Reading GTDB taxonomy for representative genomes.')
        taxonomy = Taxonomy().read(Config.TAXONOMY_FILE)

        if options.gtdbtk_classification_file:
            # add and overwrite taxonomy for genomes specified in the
            # GTDB-Tk classification file
            check_file_exists(options.gtdbtk_classification_file)

            self.logger.info('Reading GTDB-Tk classification file.')
            gtdbtk_taxonomy = Taxonomy().read(options.gtdbtk_classification_file)
            del gtdbtk_taxonomy['user_genome']
            num_reassigned = 0
            for gid, taxa in gtdbtk_taxonomy.items():
                if gid in taxonomy:
                    num_reassigned += 1
                taxonomy[gid] = taxa

            self.logger.info(f'Read GTDB-Tk classifications for {len(gtdbtk_taxonomy):,} genomes.')
            self.logger.info(f'Reassigned taxonomy for {num_reassigned:,} GTDB representative genomes.')

        if options.custom_taxonomy_file:
            # add and overwrite taxonomy for genomes specified in the
            # custom taxonomy file
            check_file_exists(options.custom_taxonomy_file)

            self.logger.info('Reading custom taxonomy file.')
            custom_taxonomy = Taxonomy().read(options.custom_taxonomy_file)
            num_reassigned = 0
            for gid, taxa in custom_taxonomy.items():
                if gid in taxonomy:
                    num_reassigned += 1
                taxonomy[gid] = taxa

            self.logger.info(f'Read custom taxonomy for {len(custom_taxonomy):,} genomes.')
            self.logger.info(f'Reassigned taxonomy for {num_reassigned:,} GTDB representative genomes.')

        if options.gtdbtk_classification_file and options.custom_taxonomy_file:
            dup_genomes = set(gtdbtk_taxonomy).intersection(custom_taxonomy)
            if len(dup_genomes) > 0:
                self.logger.error('GTDB-Tk classification and custom taxonomy '
                                  'files must not specify taxonomies for the '
                                  'same genomes.')
                self.logger.error('These files have {:,} genomes in common.'.format(len(dup_genomes)))
                self.logger.error('Example duplicate genome: {}'.format(dup_genomes.pop()))
                raise GTDBTkExit('Duplicated taxonomy information.')

        self.logger.info(f'Read taxonomy for {len(taxonomy):,} genomes.')

        return taxonomy

    def identify(self, options):
        """Identify marker genes in genomes.

        Parameters
        ----------
        options : argparse.Namespace
            The CLI arguments input by the user.
        """

        if options.genome_dir:
            check_dir_exists(options.genome_dir)

        if options.batchfile:
            check_file_exists(options.batchfile)

        make_sure_path_exists(options.out_dir)

        genomes, tln_tables = self._genomes_to_process(options.genome_dir,
                                                       options.batchfile,
                                                       options.extension)
        self.genomes_to_process = genomes

        markers = Markers(options.cpus)
        markers.identify(genomes,
                         tln_tables,
                         options.out_dir,
                         options.prefix,
                         options.force,
                         options.write_single_copy_genes)

        self.logger.info('Done.')

    def align(self, options):
        """Create MSA from marker genes.

        Parameters
        ----------
        options : argparse.Namespace
            The CLI arguments input by the user.
        """

        check_dir_exists(options.identify_dir)
        make_sure_path_exists(options.out_dir)

        markers = Markers(options.cpus, options.debug)
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
                      options.outgroup_taxon if hasattr(options, 'outgroup_taxon') else None,
                      self.genomes_to_process)

        self.logger.info('Done.')

    def infer(self, options):
        """Infer a tree from a user specified MSA.

        Parameters
        ----------
        options : argparse.Namespace
            The CLI arguments input by the user.
        """

        check_file_exists(options.msa_file)
        make_sure_path_exists(options.out_dir)

        check_dependencies(['FastTree' + ('MP' if options.cpus > 1 else '')])

        if hasattr(options, 'suffix'):
            output_tree = os.path.join(options.out_dir,
                                       PATH_MARKER_UNROOTED_TREE.format(prefix=options.prefix,
                                                                        marker=options.suffix))
            tree_log = os.path.join(options.out_dir,
                                    PATH_MARKER_TREE_LOG.format(prefix=options.prefix,
                                                                marker=options.suffix))
            fasttree_log = os.path.join(options.out_dir,
                                        PATH_MARKER_FASTTREE_LOG.format(prefix=options.prefix,
                                                                        marker=options.suffix))
        else:
            output_tree = os.path.join(options.out_dir,
                                       PATH_UNROOTED_TREE.format(prefix=options.prefix))
            tree_log = os.path.join(options.out_dir,
                                    PATH_TREE_LOG.format(prefix=options.prefix))
            fasttree_log = os.path.join(options.out_dir,
                                        PATH_FASTTREE_LOG.format(prefix=options.prefix))

        fasttree = FastTree()
        fasttree.run(output_tree, tree_log, fasttree_log, options.prot_model,
                     options.no_support, options.gamma, options.msa_file,
                     options.cpus)
        self.logger.info(f'FastTree version: {fasttree.version}')

        if hasattr(options, 'subparser_name') and options.subparser_name == 'infer':
            symlink_f(output_tree[len(options.out_dir) + 1:],
                      os.path.join(options.out_dir,
                                   os.path.basename(output_tree)))

        self.logger.info('Done.')

    def run_test(self, options):
        """Run test of classify workflow.

        Parameters
        ----------
        options : argparse.Namespace
            The CLI arguments input by the user.

        Returns
        -------
        bool
            True if the test succeeds.

        Raises
        ------
        GTDBTkTestFailure
            If the test fails.
        """

        # Use a temporary directory if none is supplied.
        if options.out_dir:
            out_dir_fh = None
            make_sure_path_exists(options.out_dir)
        else:
            out_dir_fh = tempfile.TemporaryDirectory(prefix='gtdbtk_tmp_')
            options.out_dir = out_dir_fh.name
            self.logger.info('Using a temporary directory as out_dir was not specified.')

        try:
            output_dir = os.path.join(options.out_dir, 'output')
            genome_test_dir = os.path.join(options.out_dir, 'genomes')
            if os.path.exists(genome_test_dir):
                self.logger.error(f'Test directory {genome_test_dir} already exists.')
                self.logger.error('Test must be run in a new directory.')
                sys.exit(1)

            current_path = os.path.dirname(os.path.realpath(__file__))
            input_dir = os.path.join(current_path, 'tests', 'data', 'genomes')

            shutil.copytree(input_dir, genome_test_dir)

            args = ['gtdbtk', 'classify_wf', '--genome_dir', genome_test_dir,
                    '--out_dir', output_dir, '--cpus', str(options.cpus)]
            self.logger.info('Command: {}'.format(' '.join(args)))

            # Pipe the output and write to disk.
            path_stdout = os.path.join(options.out_dir, 'test_execution.log')
            with subprocess.Popen(args, stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE, encoding='utf-8') as proc:
                with open(path_stdout, 'w') as fh_stdout:
                    bar_fmt = ' <TEST OUTPUT> '.center(22) + '{desc}'
                    with tqdm(bar_format=bar_fmt, leave=False) as p_bar:
                        while True:
                            line = proc.stdout.readline()
                            if not line:
                                break
                            fh_stdout.write(f'{line}')
                            p_bar.set_description_str(line.strip())
                proc.wait()
                exit_code = proc.returncode

            summary_fh = ClassifySummaryFileAR122(output_dir, 'gtdbtk')

            if exit_code != 0:
                self.logger.error('The test returned a non-zero exit code.')
                self.logger.error('A detailed summary of the execution log can be '
                                  'found here: {}'.format(path_stdout))
                self.logger.error('The test has failed.')
                sys.exit(1)
            if not os.path.exists(summary_fh.path):
                self.logger.error(f"{summary_fh.path} is missing.")
                self.logger.error('A detailed summary of the execution log can be '
                                  'found here: {}'.format(path_stdout))
                self.logger.error('The test has failed.')
                sys.exit(1)
        finally:
            if out_dir_fh:
                out_dir_fh.cleanup()

        self.logger.info('Test has successfully finished.')
        return True

    def classify(self, options):
        """Determine taxonomic classification of genomes.

        Parameters
        ----------
        options : argparse.Namespace
            The CLI arguments input by the user.
        """

        check_dir_exists(options.align_dir)
        make_sure_path_exists(options.out_dir)
        if options.scratch_dir:
            make_sure_path_exists(options.scratch_dir)

        genomes, _ = self._genomes_to_process(options.genome_dir,
                                              options.batchfile,
                                              options.extension)

        classify = Classify(options.cpus, options.pplacer_cpus, options.min_af)
        classify.run(genomes,
                     options.align_dir,
                     options.out_dir,
                     options.prefix,
                     options.scratch_dir,
                     options.recalculate_red,
                     options.debug,
                     options.split_tree)

        self.logger.info('Done.')

    def trim_msa(self, options):
        """ Trim an untrimmed archaea or bacterial MSA file.

        Parameters
        ----------
        options : argparse.Namespace
            The CLI arguments input by the user.
        """
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
        """Export the untrimmed archaeal or bacterial MSA file.

        Parameters
        ----------
        options : argparse.Namespace
            The CLI arguments input by the user.
        """
        misc = Misc()
        misc.export_msa(options.domain, options.output)

        self.logger.info('Done.')

    def root(self, options):
        """Root tree using outgroup.

        Parameters
        ----------
        options : argparse.Namespace
            The CLI arguments input by the user.
        """

        check_file_exists(options.input_tree)

        taxonomy = self._read_taxonomy_files(options)

        self.logger.info(f'Identifying genomes from the specified outgroup: {options.outgroup_taxon}')
        outgroup = set()
        for genome_id, taxa in taxonomy.items():
            if options.outgroup_taxon in taxa:
                outgroup.add(genome_id)

        reroot = RerootTree()
        reroot.root_with_outgroup(options.input_tree,
                                  options.output_tree,
                                  outgroup)

        self.logger.info('Done.')

    def decorate(self, options):
        """Decorate tree with GTDB taxonomy.

        Parameters
        ----------
        options : argparse.Namespace
            The CLI arguments input by the user.
        """

        check_file_exists(options.input_tree)

        taxonomy = self._read_taxonomy_files(options)

        d = Decorate()
        d.run(options.input_tree,
              taxonomy,
              options.output_tree)

        self.logger.info('Done.')

        # symlink to the decorated tree file, if not run independently
        if hasattr(options, 'suffix'):
            if options.suffix == 'bac120':
                symlink_f(PATH_BAC120_DECORATED_TREE.format(prefix=options.prefix),
                          os.path.join(options.out_dir,
                                       os.path.basename(PATH_BAC120_DECORATED_TREE.format(prefix=options.prefix))))
                symlink_f(PATH_BAC120_DECORATED_TREE.format(prefix=options.prefix) + '-table',
                          os.path.join(options.out_dir,
                                       os.path.basename(
                                           PATH_BAC120_DECORATED_TREE.format(prefix=options.prefix) + '-table')))
            elif options.suffix == 'ar122':
                symlink_f(PATH_AR122_DECORATED_TREE.format(prefix=options.prefix),
                          os.path.join(options.out_dir,
                                       os.path.basename(PATH_AR122_DECORATED_TREE.format(prefix=options.prefix))))
                symlink_f(PATH_AR122_DECORATED_TREE.format(prefix=options.prefix) + '-table',
                          os.path.join(options.out_dir,
                                       os.path.basename(
                                           PATH_AR122_DECORATED_TREE.format(prefix=options.prefix) + '-table')))
            else:
                raise GenomeMarkerSetUnknown(
                    'There was an error determining the marker set.')

    def check_install(self):
        """ Verify all GTDB-Tk data files are present.

        Raises
        ------
        ReferenceFileMalformed
            If one or more reference files are malformed.
        """
        self.logger.info("Running install verification")
        misc = Misc()
        misc.check_install()
        self.logger.info('Done.')

    def infer_ranks(self, options):
        """Establish taxonomic ranks of internal nodes using RED."""

        check_file_exists(options.input_tree)

        p = InferRanks()
        p.run(options.input_tree,
              options.ingroup_taxon,
              options.output_tree)

        self.logger.info('Done.')

    def ani_rep(self, options):
        """Calculates ANI to GTDB representative genomes.

        Parameters
        ----------
        options : argparse.Namespace
            The CLI arguments input by the user.
        """
        make_sure_path_exists(options.out_dir)

        genomes, _ = self._genomes_to_process(options.genome_dir,
                                              options.batchfile,
                                              options.extension)

        ani_rep = ANIRep(options.cpus)
        ani_rep.run(genomes, options.no_mash, options.mash_d, options.out_dir, options.prefix,
                    options.mash_k, options.mash_v, options.mash_s, options.min_af, options.mash_db)

        self.logger.info('Done.')

    def parse_options(self, options):
        """Parse user options and call the correct pipeline(s)

        Parameters
        ----------
        options : argparse.Namespace
            The CLI arguments input by the user.
        """

        # Stop processing if python 2 is being used.
        if sys.version_info.major < 3:
            raise GTDBTkExit('Python 2 is no longer supported.')

        # Correct user paths
        if hasattr(options, 'out_dir') and options.out_dir:
            options.out_dir = os.path.expanduser(options.out_dir)

        # Assert that the number of CPUs is a positive integer.
        if hasattr(options, 'cpus') and options.cpus < 1:
            self.logger.warning(
                'You cannot use less than 1 CPU, defaulting to 1.')
            options.cpus = 1

        if options.subparser_name == 'de_novo_wf':
            check_dependencies(['prodigal', 'hmmalign'])
            check_dependencies(['FastTree' + ('MP' if options.cpus > 1 else '')])

            options.write_single_copy_genes = False
            self.identify(options)

            options.identify_dir = options.out_dir
            options.skip_trimming = False
            self.align(options)

            if options.bacteria:
                options.suffix = "bac120"
            else:
                options.suffix = "ar122"

            if options.skip_gtdb_refs:
                if options.suffix == 'bac120':
                    options.msa_file = os.path.join(
                        options.out_dir, PATH_BAC120_USER_MSA.format(prefix=options.prefix))
                elif options.suffix == 'ar122':
                    options.msa_file = os.path.join(
                        options.out_dir, PATH_AR122_USER_MSA.format(prefix=options.prefix))
                else:
                    self.logger.error(
                        'There was an error determining the marker set.')
                    raise GenomeMarkerSetUnknown(
                        'Unknown marker set: {}'.format(options.suffix))
            else:
                if options.suffix == 'bac120':
                    options.msa_file = os.path.join(
                        options.out_dir, PATH_BAC120_MSA.format(prefix=options.prefix))
                elif options.suffix == 'ar122':
                    options.msa_file = os.path.join(
                        options.out_dir, PATH_AR122_MSA.format(prefix=options.prefix))
                else:
                    self.logger.error(
                        'There was an error determining the marker set.')
                    raise GenomeMarkerSetUnknown(
                        'Unknown marker set: {}'.format(options.suffix))

            self.infer(options)

            if options.suffix == 'bac120':
                options.input_tree = os.path.join(options.out_dir,
                                                  PATH_BAC120_UNROOTED_TREE.format(prefix=options.prefix))
                options.output_tree = os.path.join(options.out_dir,
                                                   PATH_BAC120_ROOTED_TREE.format(prefix=options.prefix))
            elif options.suffix == 'ar122':
                options.input_tree = os.path.join(options.out_dir,
                                                  PATH_AR122_UNROOTED_TREE.format(prefix=options.prefix))
                options.output_tree = os.path.join(options.out_dir,
                                                   PATH_AR122_ROOTED_TREE.format(prefix=options.prefix))

            self.root(options)

            if options.suffix == 'bac120':
                options.input_tree = os.path.join(options.out_dir,
                                                  PATH_BAC120_ROOTED_TREE.format(prefix=options.prefix))
                options.output_tree = os.path.join(options.out_dir,
                                                   PATH_BAC120_DECORATED_TREE.format(prefix=options.prefix))
            elif options.suffix == 'ar122':
                options.input_tree = os.path.join(options.out_dir,
                                                  PATH_AR122_ROOTED_TREE.format(prefix=options.prefix))
                options.output_tree = os.path.join(options.out_dir,
                                                   PATH_AR122_DECORATED_TREE.format(prefix=options.prefix))

            self.decorate(options)

        elif options.subparser_name == 'classify_wf':

            # TODO: Remove this block once the split_tree function is implemented.
            if hasattr(options, 'split_tree') and options.split_tree:
                self.logger.warning('The split tree option is not yet '
                                    ' supported, overriding value to False.')
            options.split_tree = False

            check_dependencies(['prodigal', 'hmmalign', 'pplacer', 'guppy',
                                'fastANI'])

            options.write_single_copy_genes = False
            self.identify(options)

            options.identify_dir = options.out_dir
            options.align_dir = options.out_dir
            options.taxa_filter = None
            options.custom_msa_filters = False
            # Added here due to the other mutex argument being include above.
            options.skip_trimming = False
            options.min_consensus = None
            options.min_perc_taxa = None
            options.skip_gtdb_refs = False
            options.cols_per_gene = None
            options.max_consensus = None
            options.rnd_seed = None
            options.skip_trimming = False
            options.scratch_dir = None
            options.recalculate_red = False

            self.align(options)

            self.classify(options)
        elif options.subparser_name == 'identify':
            self.identify(options)
        elif options.subparser_name == 'align':
            self.align(options)
        elif options.subparser_name == 'infer':
            self.infer(options)
        elif options.subparser_name == 'classify':

            # TODO: Remove this block once the split_tree function is implemented.
            if hasattr(options, 'split_tree') and options.split_tree:
                self.logger.warning('The split tree option is not yet '
                                    ' supported, overriding value to False.')
            options.split_tree = False

            if options.recalculate_red and options.split_tree:
                raise GTDBTkExit('--split_tree and --recalculate_red are mutually exclusive.')
            self.classify(options)
        elif options.subparser_name == 'root':
            self.root(options)
        elif options.subparser_name == 'decorate':
            self.decorate(options)
        elif options.subparser_name == 'infer_ranks':
            self.infer_ranks(options)
        elif options.subparser_name == 'ani_rep':
            self.ani_rep(options)
        elif options.subparser_name == 'trim_msa':
            self.trim_msa(options)
        elif options.subparser_name == 'export_msa':
            self.export_msa(options)
        elif options.subparser_name == 'test':
            check_dependencies(['prodigal', 'hmmalign', 'pplacer', 'guppy',
                                'fastANI'])
            self.run_test(options)
        elif options.subparser_name == 'check_install':
            self.check_install()
        else:
            self.logger.error('Unknown GTDB-Tk command: "' +
                              options.subparser_name + '"\n')
            sys.exit(1)

        return 0
