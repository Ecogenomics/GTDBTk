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
import argparse
import logging
import ntpath
import os
import shutil
import subprocess
import sys
import tempfile
import time
from datetime import datetime, timedelta
from typing import Dict, Tuple

from tqdm import tqdm

import gtdbtk.config.config as Config
from gtdbtk.ani_rep import ANIRep
from gtdbtk.ani_screen import ANIScreener
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
from gtdbtk.files.stage_logger import StageLoggerFile, ANIScreenStep, IdentifyStep, ClassifyStep, AlignStep, \
    InferStep, RootStep, DecorateStep
from gtdbtk.infer_ranks import InferRanks
from gtdbtk.files.batchfile import Batchfile
from gtdbtk.files.classify_summary import ClassifySummaryFileAR53, ClassifySummaryFile
from gtdbtk.markers import Markers
from gtdbtk.misc import Misc
from gtdbtk.model.enum import Domain
from gtdbtk.pipeline.export_msa import export_msa
from gtdbtk.reroot_tree import RerootTree
from gtdbtk.tools import symlink_f, get_reference_ids, confirm, assert_outgroup_taxon_valid


class OptionsParser(object):

    def __init__(self, version,output_dir=None):
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

        #Setup the Stage Logger File
        if output_dir is not None:
            base_name = ntpath.basename(sys.argv[0])
            if base_name == '__main__.py':
                prog_name = __name__.split('.')[0]
            else:
                prog_name = base_name
            #timestamp_logger.info(f'{prog_name} {" ".join(sys.argv[1:])}')
            self.stage_logger_file = StageLoggerFile(output_dir=output_dir,
                                                    version=self.version,
                                                    command_line=f'{prog_name} {" ".join(sys.argv[1:])}',
                                                    database_version = Config.VERSION_DATA,
                                                    database_path=Config.GENERIC_PATH)
            self.stage_logger = self.stage_logger_file.stage_logger

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

    def _verify_genome_id(self, genome_id: str) -> bool:
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
        GTDBTkExit
            If the genome identifier contains illegal characters.
        """
        if genome_id is None or not isinstance(genome_id, str):
            raise GTDBTkExit(f'The genome name is not a valid string: {genome_id}')
        if len(genome_id) == 0:
            raise GTDBTkExit('Genome name cannot be blank, check for input files '
                             'without a name, or empty columns in the batchfile.')
        invalid_chars = frozenset('()[],;= ')
        if any((c in invalid_chars) for c in genome_id):
            self.logger.error(f'Invalid genome ID: {genome_id}')
            self.logger.error(f'The following characters are invalid: '
                              f'{" ".join(invalid_chars)}')
            raise GTDBTkExit(f'Invalid genome ID: {genome_id}')
        return True

    @staticmethod
    def _verify_file_path(file_path: str) -> bool:
        if ' ' in file_path:
            raise GTDBTkExit(f'The genome path contains a space, this is '
                             f'unsupported by downstream applications: {file_path}')
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

        # Check that there are no illegal characters in the file path
        for file_path in genomic_files.values():
            self._verify_file_path(file_path)

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
                             f'batchfile which do not exist, see gtdbtk.warnings.log')

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
            gtdbtk_classify_file = ClassifySummaryFile(path=options.gtdbtk_classification_file)
            gtdbtk_classify_file.read()
            gtdbtk_taxonomy = gtdbtk_classify_file.get_gid_taxonomy()
            if len(gtdbtk_taxonomy) == 0:
                raise GTDBTkExit(f'No genomes found in GTDB-Tk classification file: {options.gtdbtk_classification_file}')

            num_rep_reassigned = 0
            num_usr_reassigned = 0
            for gid, taxa in gtdbtk_taxonomy.items():
                if gid in taxonomy:
                    num_rep_reassigned += 1
                else:
                    num_usr_reassigned += 1
                taxonomy[gid] = taxa

            self.logger.info(f'Read GTDB-Tk classifications for {len(gtdbtk_taxonomy):,} genomes.')
            self.logger.info(f'Reassigned taxonomy for {num_rep_reassigned:,} GTDB representative '
                             f'genomes, and {num_usr_reassigned:,} query genomes.')

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

    def identify(self, options , classified_genomes = None):
        """Identify marker genes in genomes.

        Parameters
        ----------
        options : argparse.Namespace
            The CLI arguments input by the user.
        """
        identify_step = IdentifyStep()
        identify_step.starts_at = datetime.now()
        identify_step.output_dir = options.out_dir
        identify_step.genes = options.genes
        identify_step.extension = options.extension
        identify_step.write_single_copy_genes = options.write_single_copy_genes

        if options.genome_dir:
            check_dir_exists(options.genome_dir)
            identify_step.genome_dir = options.genome_dir

        if options.batchfile:
            check_file_exists(options.batchfile)
            identify_step.batchfile = options.batchfile

        make_sure_path_exists(options.out_dir)

        genomes, tln_tables = self._genomes_to_process(options.genome_dir,
                                                       options.batchfile,
                                                       options.extension)
        #if self.genomes_to_process does not exist, create it
        if not hasattr(self, 'genomes_to_process'):
            self.genomes_to_process = genomes

        if classified_genomes is not None:
            for marker_set_id in ['bac120', 'ar53']:
                if marker_set_id in classified_genomes:
                    for classified_genome in classified_genomes.get(marker_set_id).keys():
                        genomes.pop(classified_genome, None)

        markers = Markers(options.cpus)
        reports = markers.identify(genomes,
                         tln_tables,
                         options.out_dir,
                         options.prefix,
                         options.force,
                         options.genes,
                         options.write_single_copy_genes)

        identify_step.output_files = reports

        identify_step.ends_at = datetime.now()
        duration = identify_step.ends_at - identify_step.starts_at
        #we round the duration to the nearest second
        identify_step.duration = str(duration - timedelta(microseconds=duration.microseconds))
        identify_step.status = 'completed'

        self.stage_logger.steps.append(identify_step)
        self.logger.info('Done.')

    def align(self, options):
        """Create MSA from marker genes.

        Parameters
        ----------
        options : argparse.Namespace
            The CLI arguments input by the user.
        """

        align_step = AlignStep()
        align_step.starts_at = datetime.now()
        align_step.output_dir = options.out_dir
        align_step.skip_gtdb_refs = options.skip_gtdb_refs
        align_step.taxa_filter = options.taxa_filter
        align_step.min_perc_aa = options.min_perc_aa
        align_step.custom_msa_filters = options.custom_msa_filters
        align_step.skip_trimming = options.skip_trimming
        align_step.rnd_seed = options.rnd_seed
        align_step.cols_per_gene = options.cols_per_gene
        align_step.min_consensus = options.min_consensus
        align_step.max_consensus = options.max_consensus
        align_step.min_perc_taxa = options.min_perc_taxa
        align_step.outgroup_taxon = options.outgroup_taxon if hasattr(options, 'outgroup_taxon') else None

        check_dir_exists(options.identify_dir)
        align_step.identify_dir = options.identify_dir
        make_sure_path_exists(options.out_dir)

        markers = Markers(options.cpus, options.debug)
        reports = markers.align(options.identify_dir,
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

        align_step.ends_at = datetime.now()
        duration = align_step.ends_at - align_step.starts_at
        #we round the duration to the nearest second
        align_step.duration = str(duration - timedelta(microseconds=duration.microseconds))
        align_step.output_files = reports
        align_step.status = 'completed'

        self.stage_logger.steps.append(align_step)

        self.logger.info('Done.')


    def infer(self, options):
        """Infer a tree from a user specified MSA.

        Parameters
        ----------
        options : argparse.Namespace
            The CLI arguments input by the user.
        """

        infer_step = InferStep()
        infer_step.starts_at = datetime.now()
        infer_step.output_dir = options.out_dir
        infer_step.msa_file = options.msa_file
        infer_step.prot_model = options.prot_model
        infer_step.no_support = options.no_support
        infer_step.gamma = options.gamma


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
            symlink_f(output_tree[len(options.out_dir.rstrip('/')) + 1:],
                      os.path.join(options.out_dir,
                                   os.path.basename(output_tree)))

        infer_step.ends_at = datetime.now()
        duration = infer_step.ends_at - infer_step.starts_at
        #we round the duration to the nearest second
        infer_step.duration = str(duration - timedelta(microseconds=duration.microseconds))
        infer_step.status = 'completed'

        self.stage_logger.steps.append(infer_step)

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
                    '--out_dir', output_dir, '--cpus', str(options.cpus), '-f']
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

            summary_fh = ClassifySummaryFileAR53(output_dir, 'gtdbtk')

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

        classify_step = ClassifyStep()
        classify_step.starts_at = datetime.now()
        classify_step.output_dir = options.out_dir
        classify_step.debug_option = options.debug
        classify_step.full_tree = options.full_tree
        classify_step.skip_ani_screen = options.skip_ani_screen
        classify_step.no_mash = options.no_mash
        classify_step.mash_k = options.mash_k
        classify_step.mash_v = options.mash_v
        classify_step.mash_s = options.mash_s
        classify_step.mash_db = options.mash_db
        classify_step.mash_max_dist = options.mash_max_distance

        ani_summary_files = {}
        if self.stage_logger_file.stage_logger.has_stage(ANIScreenStep):
            previous_ani_step = self.stage_logger_file.stage_logger.get_stage(ANIScreenStep)
            ani_summary_files = previous_ani_step.output_files


        check_dir_exists(options.align_dir)
        classify_step.align_dir = options.align_dir
        make_sure_path_exists(options.out_dir)
        if options.scratch_dir:
            make_sure_path_exists(options.scratch_dir)
            classify_step.scratch_dir = options.scratch_dir

        if options.genome_dir:
            check_dir_exists(options.genome_dir)
            classify_step.genome_dir = options.genome_dir

        if options.batchfile:
            check_file_exists(options.batchfile)
            classify_step.batchfile = options.batchfile

        genomes, _ = self._genomes_to_process(options.genome_dir,
                                              options.batchfile,
                                              options.extension)

        classify = Classify(options.cpus, options.pplacer_cpus, options.min_af)
        reports = classify.run(genomes=genomes,
                     align_dir=options.align_dir,
                     out_dir=options.out_dir,
                     prefix=options.prefix,
                     scratch_dir=options.scratch_dir,
                     debugopt=options.debug,
                     fulltreeopt=options.full_tree,
                     skip_ani_screen=options.skip_ani_screen,
                     no_mash=options.no_mash,
                     mash_k=options.mash_k,
                     mash_v=options.mash_v,
                     mash_s=options.mash_s,
                     mash_db=options.mash_db,
                     mash_max_dist=options.mash_max_distance,
                     ani_summary_files=ani_summary_files
                     )

        classify_step.ends_at = datetime.now()
        duration = classify_step.ends_at - classify_step.starts_at
        classify_step.output_files = reports
        #we round the duration to the nearest second
        classify_step.duration = str(duration - timedelta(microseconds=duration.microseconds))
        classify_step.status = 'completed'

        self.stage_logger.steps.append(classify_step)

        self.logger.info('Note that Tk classification mode is insufficient for publication of new taxonomic '
                         'designations. New designations should be based on one or more de novo trees, an '
                         'example of which can be produced by Tk in de novo mode.')

        self.logger.info('Done.')

    def ani_screen(self, options ):
        """Run a mash/FastANI screen of all user genomes
        against the reference genomes.

        Parameters
        ----------
        options : argparse.Namespace
            The CLI arguments input by the user.
        """

        ani_step = ANIScreenStep()
        ani_step.starts_at = datetime.now()
        ani_step.output_dir = options.out_dir
        ani_step.mash_db = options.mash_db
        ani_step.mash_k = options.mash_k
        ani_step.mash_v = options.mash_v
        ani_step.mash_s = options.mash_s
        ani_step.mash_max_dist = options.mash_max_distance

        if options.genome_dir:
            check_dir_exists(options.genome_dir)
            ani_step.genome_dir = options.genome_dir

        if options.batchfile:
            check_file_exists(options.batchfile)
            ani_step.batchfile = options.batchfile

        make_sure_path_exists(options.out_dir)

        genomes, tln_tables = self._genomes_to_process(options.genome_dir,
                                                       options.batchfile,
                                                       options.extension)
        self.genomes_to_process = genomes

        make_sure_path_exists(options.out_dir)

        genomes, _ = self._genomes_to_process(options.genome_dir,
                                              options.batchfile,
                                              options.extension)

        aniscreener = ANIScreener(options.cpus)
        classified_genomes,reports = aniscreener.run_aniscreen(
            genomes=genomes,
            no_mash=options.no_mash,
            out_dir=options.out_dir,
            prefix=options.prefix,
            mash_k=options.mash_k,
            mash_v=options.mash_v,
            mash_s=options.mash_s,
            mash_max_dist=options.mash_max_distance,
            mash_db=options.mash_db)


        self.logger.info('Done.')

        ani_step.ends_at = datetime.now()
        duration = ani_step.ends_at - ani_step.starts_at
        #we round the duration to the nearest second
        ani_step.duration = str(duration - timedelta(microseconds=duration.microseconds))
        ani_step.status = 'completed'
        ani_step.output_files=reports

        self.stage_logger.steps.append(ani_step)

        return classified_genomes

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
        # Get and validate the domain
        user_domain = getattr(options, 'domain', None)
        if user_domain == 'arc':
            domain = Domain.ARCHAEA
        elif user_domain == 'bac':
            domain = Domain.BACTERIA
        else:
            raise GTDBTkExit(f'Unknown domain: {options.domain}')

        # Get and validate the path
        user_output = getattr(options, 'output', None)
        if not user_output:
            raise GTDBTkExit(f'You must specify a valid path: "{user_output}"')

        # Export the MSA
        export_msa(domain=domain, output_file=user_output)
        self.logger.info('Done.')

    def root(self, options: argparse.Namespace):
        """Root tree using outgroup.

        Parameters
        ----------
        options : argparse.Namespace
            The CLI arguments input by the user.
        """

        root_step = RootStep()
        root_step.starts_at = datetime.now()


        check_file_exists(options.input_tree)
        assert_outgroup_taxon_valid(options.outgroup_taxon)

        root_step.outgroup_taxon = options.outgroup_taxon

        root_step.input_tree = options.input_tree
        # set gtdbtk_classification_file exists in options we can use it
        root_step.gtdbtk_classification_file=options.gtdbtk_classification_file if hasattr(options, 'gtdbtk_classification_file') else None
        root_step.custom_taxonomy_file=options.custom_taxonomy_file if hasattr(options, 'custom_taxonomy_file') else None

        taxonomy = self._read_taxonomy_files(options)

        self.logger.info(f'Identifying genomes from the specified outgroup: {options.outgroup_taxon}')
        outgroup = set()
        for genome_id, taxa in taxonomy.items():
            if options.outgroup_taxon in taxa:
                outgroup.add(genome_id)

        reroot = RerootTree()
        reports = reroot.root_with_outgroup(options.input_tree,
                                  options.output_tree,
                                  outgroup)

        root_step.ends_at = datetime.now()
        duration = root_step.ends_at - root_step.starts_at
        #we round the duration to the nearest second
        root_step.duration = str(duration - timedelta(microseconds=duration.microseconds))
        root_step.status = 'completed'

        # if self.stage_logger exists, we add the step to the stage_logger
        if hasattr(self, 'stage_logger'):
            self.stage_logger.steps.append(root_step)

        self.logger.info('Done.')

    def decorate(self, options):
        """Decorate tree with GTDB taxonomy.

        Parameters
        ----------
        options : argparse.Namespace
            The CLI arguments input by the user.
        """

        decorate_step = DecorateStep()
        decorate_step.starts_at = datetime.now()

        check_file_exists(options.input_tree)

        taxonomy = self._read_taxonomy_files(options)

        decorate_step.input_tree = options.input_tree
        # set gtdbtk_classification_file exists in options we can use it
        decorate_step.gtdbtk_classification_file=options.gtdbtk_classification_file if hasattr(options, 'gtdbtk_classification_file') else None
        decorate_step.custom_taxonomy_file=options.custom_taxonomy_file if hasattr(options, 'custom_taxonomy_file') else None
        decorate_step.suffix=options.suffix if hasattr(options, 'suffix') else None

        d = Decorate()
        reports = d.run(options.input_tree,
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
            elif options.suffix == 'ar53':
                symlink_f(PATH_AR53_DECORATED_TREE.format(prefix=options.prefix),
                          os.path.join(options.out_dir,
                                       os.path.basename(PATH_AR53_DECORATED_TREE.format(prefix=options.prefix))))
                symlink_f(PATH_AR53_DECORATED_TREE.format(prefix=options.prefix) + '-table',
                          os.path.join(options.out_dir,
                                       os.path.basename(
                                           PATH_AR53_DECORATED_TREE.format(prefix=options.prefix) + '-table')))
            else:
                raise GenomeMarkerSetUnknown(
                    'There was an error determining the marker set.')

        decorate_step.ends_at = datetime.now()
        duration = decorate_step.ends_at - decorate_step.starts_at
        #we round the duration to the nearest second
        decorate_step.duration = str(duration - timedelta(microseconds=duration.microseconds))
        decorate_step.status = 'completed'
        decorate_step.output_files=reports

        if hasattr(self, 'stage_logger'):
            self.stage_logger.steps.append(decorate_step)

        self.logger.info('Done.')

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

    def convert_to_itol(self, options):
        """Convert Tree to iTOL format.

        Parameters
        ----------
        options : argparse.Namespace
            The CLI arguments input by the user.
        """
        check_file_exists(options.input_tree)

        r = Misc()
        r.convert_to_itol(options.input_tree, options.output_tree)
        self.logger.info('Done.')

    def remove_labels(self, options):
        """Remove labels from tree.

        Parameters
        ----------
        options : argparse.Namespace
            The CLI arguments input by the user.
        """

        check_file_exists(options.input_tree)

        r = Misc()
        r.remove_labels(options.input_tree, options.output_tree)
        self.logger.info('Done.')

    def remove_intermediate_files(self,out_dir,workflow_name):
        """Remove intermediate files from the output directory.
        Parameters
        ----------
            out_dir : str
                The output directory.
        """

        misc = Misc()
        misc.remove_intermediate_files(out_dir,workflow_name)
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
            assert_outgroup_taxon_valid(options.outgroup_taxon)

            if options.skip_gtdb_refs and options.custom_taxonomy_file is None:
                raise GTDBTkExit("When running de_novo_wf, The '--skip_gtdb_refs' flag requires"
                                 "'--custom_taxonomy_file' to be included to the command line.")

            if options.write_single_copy_genes and not options.keep_intermediates:
                self.logger.warning('--write_single_copy_genes flag is set to True,'
                                    ' but --keep_intermediates is set to False. '
                                    'The intermediate folder containing the single copy genes will be removed.')
                if not confirm('Do you want to proceed?'):
                    self.logger.info('Exiting workflow.')
                    sys.exit(0)

            self.identify(options)

            options.identify_dir = options.out_dir
            options.skip_trimming = False
            self.align(options)

            if options.bacteria:
                options.suffix = "bac120"
            else:
                options.suffix = "ar53"

            if options.skip_gtdb_refs:
                if options.suffix == 'bac120':
                    if os.path.isfile(os.path.join(options.out_dir,
                                                   PATH_BAC120_USER_MSA.format(prefix=options.prefix))):
                        options.msa_file = os.path.join(options.out_dir,
                                                     PATH_BAC120_USER_MSA.format(prefix=options.prefix))
                    else:
                        options.msa_file = os.path.join(options.out_dir,
                                                     PATH_BAC120_USER_MSA.format(prefix=options.prefix) + '.gz')

                elif options.suffix == 'ar53':
                    if os.path.isfile(os.path.join(options.out_dir,
                                                   PATH_AR53_USER_MSA.format(prefix=options.prefix))):
                        options.msa_file = os.path.join(options.out_dir,
                                                     PATH_AR53_USER_MSA.format(prefix=options.prefix))
                    else:
                        options.msa_file = os.path.join(options.out_dir,
                                                     PATH_AR53_USER_MSA.format(prefix=options.prefix) + '.gz')
                else:
                    self.logger.error(
                        'There was an error determining the marker set.')
                    raise GenomeMarkerSetUnknown(
                        'Unknown marker set: {}'.format(options.suffix))
            else:
                if options.suffix == 'bac120':
                    if os.path.isfile(os.path.join(
                        options.out_dir, PATH_BAC120_MSA.format(prefix=options.prefix))):
                        options.msa_file = os.path.join(
                            options.out_dir, PATH_BAC120_MSA.format(prefix=options.prefix))
                    else:
                        options.msa_file = os.path.join(
                            options.out_dir, PATH_BAC120_MSA.format(prefix=options.prefix) + '.gz')
                elif options.suffix == 'ar53':
                    if os.path.isfile(os.path.join(
                        options.out_dir, PATH_AR53_MSA.format(prefix=options.prefix))):
                        options.msa_file = os.path.join(
                            options.out_dir, PATH_AR53_MSA.format(prefix=options.prefix))
                    else:
                        options.msa_file = os.path.join(
                            options.out_dir, PATH_AR53_MSA.format(prefix=options.prefix) + '.gz')
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
            elif options.suffix == 'ar53':
                options.input_tree = os.path.join(options.out_dir,
                                                  PATH_AR53_UNROOTED_TREE.format(prefix=options.prefix))
                options.output_tree = os.path.join(options.out_dir,
                                                   PATH_AR53_ROOTED_TREE.format(prefix=options.prefix))

            self.root(options)

            if options.suffix == 'bac120':
                options.input_tree = os.path.join(options.out_dir,
                                                  PATH_BAC120_ROOTED_TREE.format(prefix=options.prefix))
                options.output_tree = os.path.join(options.out_dir,
                                                   PATH_BAC120_DECORATED_TREE.format(prefix=options.prefix))
            elif options.suffix == 'ar53':
                options.input_tree = os.path.join(options.out_dir,
                                                  PATH_AR53_ROOTED_TREE.format(prefix=options.prefix))
                options.output_tree = os.path.join(options.out_dir,
                                                   PATH_AR53_DECORATED_TREE.format(prefix=options.prefix))

            self.decorate(options)

            if not options.keep_intermediates:
                self.remove_intermediate_files(options.out_dir,'de_novo_wf')

        elif options.subparser_name == 'classify_wf':

            check_dependencies(['prodigal', 'hmmalign', 'pplacer', 'guppy',
                                'fastANI'])

            if options.write_single_copy_genes and not options.keep_intermediates:
                self.logger.warning('--write_single_copy_genes flag is set to True,'
                                    ' but --keep_intermediates is set to False. '
                                    'The intermediate folder containing the single copy genes will be removed.')
                if not confirm('Do you want to proceed?'):
                    self.logger.info('Exiting workflow.')
                    sys.exit(0)

            # if options.skip_aniscreen is false,
            # we need to make sure the options.mash_db is selected too to point to the folder
            # where the sketch file is.
            if not options.skip_ani_screen and not options.no_mash and not options.mash_db:
                self.logger.error('You must specify a path to the mash database with --mash_db')
                sys.exit(1)

            #options.write_single_copy_genes = False

            classified_genomes = None
            #We ned to check if the ani screen has already be ran, if so, we need to skip it.
            #Check is the gtdbtk.json file exists in the output folder
            if os.path.isfile(self.stage_logger_file.path) and 1==2:
                #If the file exists, we need to check if the ani_screen step has been ran
                #If the ani_screen step has been ran, we need to skip it.
                self.stage_logger_file.read()
                stage_logger = self.stage_logger_file.stage_logger
                if stage_logger.has_stage(ANIScreenStep):
                    # we get the genomes already classified by the ani_screen step
                    previous_ani_step = stage_logger.get_stage(ANIScreenStep)
                    if previous_ani_step.is_complete():
                        self.logger.warning('The ani_screen step has already been completed, we load existing results.')
                        ani_summary_files = previous_ani_step.output_files
                        classify_method = Classify()
                        classified_genomes=classify_method.load_fastani_results_pre_pplacer(ani_summary_files)
                        classified_genomes = classify_method.convert_rows_to_dict(classified_genomes)
                        len_mash_classified_bac120 = len(classified_genomes['bac120']) \
                            if 'bac120' in classified_genomes else 0

                        len_mash_classified_ar53 = len(classified_genomes['ar53']) \
                            if 'ar53' in classified_genomes else 0

                        self.logger.info(f'{len_mash_classified_ar53 + len_mash_classified_bac120} genome(s) have '
                                         f'been classified using the ANI pre-screening step.')

                        options.skip_ani_screen = True

            #Before identify step, we run the ani_screen step, these genomes classify with this step will not continue
            # in the identify step.

            if not options.skip_ani_screen:
                classified_genomes = self.ani_screen(options)

            self.identify(options,classified_genomes)

            options.identify_dir = options.out_dir
            options.align_dir = options.out_dir
            options.taxa_filter = None
            options.custom_msa_filters = False
            # Added here due to the other mutex argument being included above.
            options.skip_trimming = False
            options.min_consensus = None
            options.min_perc_taxa = None
            options.skip_gtdb_refs = False
            options.cols_per_gene = None
            options.max_consensus = None
            options.rnd_seed = None
            options.skip_trimming = False

            self.align(options)

            # because we run ani_screen before the identify step, we do not need to rerun the
            #ani step again, so we set the skip_aniscreen to True
            options.skip_ani_screen = True

            self.classify(options)
            if not options.keep_intermediates:
                self.remove_intermediate_files(options.out_dir,'classify_wf')


        elif options.subparser_name == 'identify':
            self.identify(options)
        elif options.subparser_name == 'align':
            self.align(options)
        elif options.subparser_name == 'infer':
            self.infer(options)
        elif options.subparser_name == 'classify':
            # if options.skip_ani_screen is not selected,
            # we need to make sure the options.mash_db is selected too to point to the folder
            # where the sketch file is.
            if not options.skip_ani_screen and not options.no_mash and not options.mash_db:
                print(options.skip_ani_screen, options.no_mash, options.mash_db)
                self.logger.error('You must specify a path to the mash database with --mash_db')
            self.classify(options)
            self.stage_logger_file.write()

        elif options.subparser_name == 'root':
            self.root(options)
        elif options.subparser_name == 'decorate':
            self.decorate(options)
        elif options.subparser_name == 'infer_ranks':
            self.infer_ranks(options)
        elif options.subparser_name == 'ani_rep':
            self.ani_rep(options)
        elif options.subparser_name == 'remove_labels':
            self.remove_labels(options)
        elif options.subparser_name == 'convert_to_itol':
            self.convert_to_itol(options)
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

        if hasattr(self,'stage_logger' ) and self.stage_logger.steps:
            self.stage_logger_file.write()

        return 0
