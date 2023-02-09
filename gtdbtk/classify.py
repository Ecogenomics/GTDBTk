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
import random
import shutil
import sys
from collections import defaultdict
from operator import itemgetter

import dendropy

from numpy import median as np_median

import gtdbtk.config.config as Config
from gtdbtk.ani_rep import ANIRep, ANISummaryFile
from gtdbtk.biolib_lite.common import make_sure_path_exists,canonical_gid
from gtdbtk.biolib_lite.execute import check_dependencies
from gtdbtk.biolib_lite.newick import parse_label
from gtdbtk.biolib_lite.seq_io import read_seq, read_fasta
from gtdbtk.biolib_lite.taxonomy import Taxonomy
from gtdbtk.config.output import *
from gtdbtk.exceptions import GenomeMarkerSetUnknown, GTDBTkExit
from gtdbtk.external.fastani import FastANI
from gtdbtk.external.pplacer import Pplacer

from gtdbtk.files.classify_summary import ClassifySummaryFileAR53, ClassifySummaryFileBAC120, ClassifySummaryFileRow
from gtdbtk.files.marker.copy_number import CopyNumberFileAR53, CopyNumberFileBAC120
from gtdbtk.files.pplacer_classification import PplacerClassifyFileBAC120, PplacerClassifyFileAR53, \
    PplacerLowClassifyFileBAC120
from gtdbtk.files.prodigal.tln_table_summary import TlnTableSummaryFile
from gtdbtk.files.red_dict import REDDictFileAR53, REDDictFileBAC120
from gtdbtk.files.gtdb_radii import GTDBRadiiFile
from gtdbtk.files.missing_genomes import DisappearingGenomesFileAR53, DisappearingGenomesFileBAC120
from gtdbtk.files.tree_mapping import GenomeMappingFile, GenomeMappingFileRow
from gtdbtk.markers import Markers
from gtdbtk.relative_distance import RelativeDistance
from gtdbtk.split import Split
from gtdbtk.tools import add_ncbi_prefix, symlink_f, get_memory_gb, get_reference_ids, TreeTraversal, \
    calculate_patristic_distance, tqdm_log, standardise_taxonomy, limit_rank, aa_percent_msa


sys.setrecursionlimit(15000)


class Classify(object):
    """Determine taxonomic classification of genomes by ML placement."""

    def __init__(self, cpus=1, pplacer_cpus=None, af_threshold=None):
        """Initialize."""

        check_dependencies(['pplacer', 'guppy', 'fastANI'])

        self.taxonomy_file = Config.TAXONOMY_FILE
        self.af_threshold = af_threshold if af_threshold else Config.AF_THRESHOLD
        self.gtdb_taxonomy = Taxonomy().read(self.taxonomy_file)

        self.order_rank = ["d__", "p__", "c__", "o__", 'f__', 'g__', 's__']

        self.DOMAIN_IDX = 0
        self.PHYLUM_IDX = 1
        self.CLASS_IDX = 2
        self.ORDER_IDX = 3
        self.FAMILY_IDX = 4

        self.logger = logging.getLogger('timestamp')

        self.cpus = max(cpus, 1)
        self.pplacer_cpus = max(pplacer_cpus if pplacer_cpus else cpus, 1)

        self.gtdb_radii = GTDBRadiiFile()

        # pplacer seems to hang if more than 64 threads are given, warn the user
        # otherwise let them if they specify pplacer_cpus.
        if pplacer_cpus and pplacer_cpus > 64:
            self.logger.warning('pplacer is known to hang if >64 CPUs are '
                                'used. You have overridden this using: '
                                '--pplacer_cpus')
        elif not pplacer_cpus and cpus > 64:
            self.logger.warning('Setting pplacer CPUs to 64, as pplacer is '
                                'known to hang if >64 are used. You can '
                                'override this using: --pplacer_cpus')
            self.pplacer_cpus = 64

        self.species_radius = self.parse_radius_file()
        self.reference_ids = get_reference_ids()

        # rank_of_interest determine the rank in the tree_mapping file for
        # lower classification
        self.rank_of_interest = "c__"

    @staticmethod
    def parse_radius_file():
        results = {}
        with open(Config.RADII_FILE) as f:
            for line in f:
                infos = line.strip().split('\t')
                gid = infos[1]
                if infos[1].startswith('GB_') or infos[1].startswith('RS_'):
                    gid = gid[3:]
                results[gid] = float(infos[2])
        return results

    def parse_leaf_to_dir_path(self, genome_id):
        """ Convert a genome id to a path.
         i.e. GCA_123456789.0 would be converted to GCA/123/456/789/

        Parameters
        ----------
        genome_id: str
            NCBI genome id GCF/GCA_xxxxxxxxx
        :return: str
            path to the genome id path
        """
        try:
            genome_path = '/'.join([genome_id[0:3], genome_id[4:7],
                                    genome_id[7:10], genome_id[10:13]])
            return genome_path
        except IndexError:
            logger = logging.getLogger('timestamp')
            logger.error('Specified path could not be created for reference genome: ' + genome_id)
            raise GTDBTkExit('Specified path could not be created for reference genome: ' + genome_id)

    def place_genomes(self,
                      user_msa_file,
                      marker_set_id,
                      out_dir,
                      prefix,
                      scratch_dir=None,
                      levelopt=None,
                      tree_iter=None,
                      number_low_trees=None,
                      idx_tree=None):
        """Place genomes into reference tree using pplacer."""

        # Warn if the memory is insufficient
        mem_warning = 'pplacer requires ~{req_gb} GB of RAM to fully load the ' \
                      '{domain} tree into memory. However, {cur_gb:,} GB was ' \
                      'detected. This may affect pplacer performance, or fail' \
                      ' if there is insufficient swap space.'
        mem_gb = get_memory_gb()
        if mem_gb is not None:
            mem_total = mem_gb['MemTotal']
            if marker_set_id == 'bac120' and levelopt is None and mem_total < Config.PPLACER_MIN_RAM_BAC_FULL:
                self.logger.warning(mem_warning.format(req_gb=Config.PPLACER_MIN_RAM_BAC_FULL,
                                                       domain='bacterial',
                                                       cur_gb=mem_total))
            if marker_set_id == 'bac120' and levelopt == 'high' and mem_total < Config.PPLACER_MIN_RAM_BAC_SPLIT:
                self.logger.warning(mem_warning.format(req_gb=Config.PPLACER_MIN_RAM_BAC_SPLIT,
                                                       domain='bacterial',
                                                       cur_gb=mem_total))
            elif marker_set_id == 'ar53' and mem_total < Config.PPLACER_MIN_RAM_ARC:
                self.logger.warning(mem_warning.format(req_gb=Config.PPLACER_MIN_RAM_ARC,
                                                       domain='archaeal',
                                                       cur_gb=mem_total))

        # rename user MSA file for compatibility with pplacer
        if not user_msa_file.endswith('.fasta') and not user_msa_file.endswith('.gz'):
            if marker_set_id == 'bac120':
                t = PATH_BAC120_USER_MSA.format(prefix=prefix)
            elif marker_set_id == 'ar53':
                t = PATH_AR53_USER_MSA.format(prefix=prefix)
            else:
                raise GenomeMarkerSetUnknown('There was an error determining the marker set.')

            shutil.copyfile(user_msa_file, t)
            user_msa_file = t

        # run pplacer to place bins in reference genome tree
        num_genomes = sum([1 for _seq_id, _seq in read_seq(user_msa_file)])

        # check if a scratch file is to be created
        pplacer_mmap_file = None
        if scratch_dir:
            self.logger.info('Using a scratch file for pplacer allocations. '
                             'This decreases memory usage and performance.')
            pplacer_mmap_file = os.path.join(
                scratch_dir, prefix + ".pplacer.scratch")
            make_sure_path_exists(scratch_dir)

        # get path to pplacer reference package
        pplacer_ref_pkg = None
        if marker_set_id == 'bac120':
            if levelopt is None:
                self.logger.log(Config.LOG_TASK,
                                f'Placing {num_genomes:,} bacterial genomes '
                                f'into reference tree with pplacer using '
                                f'{self.pplacer_cpus} CPUs (be patient).')
                pplacer_ref_pkg = os.path.join(Config.PPLACER_DIR,
                                               Config.PPLACER_BAC120_REF_PKG)

            elif levelopt == 'high':
                self.logger.log(Config.LOG_TASK,
                                f'Placing {num_genomes:,} bacterial genomes '
                                f'into backbone reference tree with pplacer using '
                                f'{self.pplacer_cpus} CPUs (be patient).')
                pplacer_ref_pkg = os.path.join(Config.BACKBONE_PPLACER_DIR,
                                               Config.BACKBONE_PPLACER_REF_PKG)
            elif levelopt == 'low':
                self.logger.log(Config.LOG_TASK,
                                f'Placing {num_genomes:,} bacterial genomes '
                                f'into class-level reference tree {tree_iter} ({idx_tree}/{number_low_trees}) with '
                                f'pplacer using {self.pplacer_cpus} CPUs '
                                f'(be patient).')
                pplacer_ref_pkg = os.path.join(Config.CLASS_LEVEL_PPLACER_DIR,
                                               Config.CLASS_LEVEL_PPLACER_REF_PKG.format(iter=tree_iter))
        elif marker_set_id == 'ar53':
            self.logger.log(Config.LOG_TASK,
                            f'Placing {num_genomes:,} archaeal genomes into '
                            f'reference tree with pplacer using '
                            f'{self.pplacer_cpus} CPUs (be patient).')
            pplacer_ref_pkg = os.path.join(Config.PPLACER_DIR,
                                           Config.PPLACER_AR53_REF_PKG)
        else:
            raise GenomeMarkerSetUnknown(f'Unknown marker set: {marker_set_id}')

        # create pplacer output directory
        pplacer_out_dir = os.path.join(out_dir, DIR_PPLACER)
        if not os.path.exists(pplacer_out_dir):
            os.makedirs(pplacer_out_dir)

        # run pplacer
        pplacer_json_out, pplacer_out = None, None
        if marker_set_id == 'bac120':
            if levelopt is None:
                pplacer_out = os.path.join(out_dir, PATH_BAC120_PPLACER_OUT)
                pplacer_json_out = os.path.join(
                    out_dir, PATH_BAC120_PPLACER_JSON)
            elif levelopt == 'high':
                pplacer_out = os.path.join(
                    out_dir, PATH_BACKBONE_BAC120_PPLACER_OUT)
                pplacer_json_out = os.path.join(
                    out_dir, PATH_BACKBONE_BAC120_PPLACER_JSON)
            elif levelopt == 'low':
                pplacer_out = os.path.join(
                    out_dir, PATH_CLASS_LEVEL_BAC120_PPLACER_OUT.format(iter=tree_iter))
                pplacer_json_out = os.path.join(
                    out_dir, PATH_CLASS_LEVEL_BAC120_PPLACER_JSON.format(iter=tree_iter))
        elif marker_set_id == 'ar53':
            pplacer_out = os.path.join(out_dir, PATH_AR53_PPLACER_OUT)
            pplacer_json_out = os.path.join(out_dir, PATH_AR53_PPLACER_JSON)
        else:
            self.logger.error('There was an error determining the marker set.')
            raise GenomeMarkerSetUnknown

        pplacer = Pplacer()
        if levelopt is None or levelopt == 'high':
            self.logger.info(f'pplacer version: {pplacer.version}')
        # #DEBUG: Skip pplacer
        run_pplacer = False
        if run_pplacer:
            pplacer.run(self.pplacer_cpus, 'wag', pplacer_ref_pkg, pplacer_json_out,
                        user_msa_file, pplacer_out, pplacer_mmap_file)

        # extract tree
        tree_file = None
        if marker_set_id == 'bac120':
            if levelopt is None:
                tree_file = os.path.join(
                    out_dir, PATH_BAC120_TREE_FILE.format(prefix=prefix))
            elif levelopt == 'high':
                tree_file = os.path.join(
                    out_dir, PATH_BACKBONE_BAC120_TREE_FILE.format(prefix=prefix))
            elif levelopt == 'low':
                tree_file = os.path.join(
                    out_dir, PATH_CLASS_LEVEL_BAC120_TREE_FILE.format(prefix=prefix, iter=tree_iter))
        elif marker_set_id == 'ar53':
            tree_file = os.path.join(
                out_dir, PATH_AR53_TREE_FILE.format(prefix=prefix))
        else:
            self.logger.error('There was an error determining the marker set.')
            raise GenomeMarkerSetUnknown

        pplacer.tog(pplacer_json_out, tree_file)
        return tree_file

    def _parse_red_dict(self, red_dist_dict):
        results = {}
        for k, v in red_dist_dict.items():
            if k in ['d__', 'domain']:
                results['d__'] = v
            elif k in ['p__', 'phylum']:
                results['p__'] = v
            elif k in ['c__', 'class']:
                results['c__'] = v
            elif k in ['o__', 'order']:
                results['o__'] = v
            elif k in ['f__', 'family']:
                results['f__'] = v
            elif k in ['g__', 'genus']:
                results['g__'] = v
            elif k in ['s__', 'species']:
                results['s__'] = v
        return results

    def parser_marker_summary_file(self, marker_summary_fh):
        results = dict()
        for gid, marker_dict in marker_summary_fh.genomes.items():
            multi_hits_percent = (100 * len(marker_dict['mul'])) / \
                                 len(marker_summary_fh.marker_names)
            if multi_hits_percent >= Config.DEFAULT_MULTIHIT_THRESHOLD:
                results[gid] = round(multi_hits_percent, 1)
        return results

    def run(self,
            genomes,
            align_dir,
            out_dir,
            prefix,
            scratch_dir=None,
            debugopt=False,
            fulltreeopt=False,
            skip_ani_screen=False,
            no_mash=False,
            mash_k=Config.MASH_K_VALUE,
            mash_v=Config.MASH_V_VALUE,
            mash_s=Config.MASH_S_VALUE,
            mash_max_dist=Config.MASH_MAX_DISTANCE,
            mash_db=None,
            ani_summary_files=None):
        """Classify genomes based on position in reference tree."""

        _bac_gids, _ar_gids, bac_ar_diff = Markers().genome_domain(align_dir, prefix)


        # If prescreen is set to True, then we will first run all genomes against a mash database
        # of all genomes in the reference package. The next step will be to classify those genomes with
        # FastANI.
        # All genomes classified with FastANI will be removed from the input genomes list for the
        # rest of the pipeline.
        mash_classified_user_genomes = {}
        if not skip_ani_screen:
            if not no_mash:
                # if mash_db finishes with a backslash, it should be considered a directory
                if mash_db.endswith('/'):
                    make_sure_path_exists(mash_db)
                if os.path.isdir(mash_db):
                    mash_db = os.path.join(mash_db, Config.MASH_SKETCH_FILE)

            # we set mash_d == mash_max_dist to avoid user to run mash with impossible values
            mash_d = mash_max_dist

            ani_rep = ANIRep(self.cpus)
            # we store all the mash information in the classify directory
            fastani_results = ani_rep.run_mash_fastani(genomes, no_mash, mash_d, os.path.join(out_dir, DIR_ANISCREEN), prefix, mash_k, mash_v, mash_s,mash_max_dist, mash_db )

            mash_classified_user_genomes = self._sort_fastani_results_pre_pplacer(
                fastani_results,bac_ar_diff)

            len_mash_classified_bac120 = len(mash_classified_user_genomes['bac120']) \
                if 'bac120' in mash_classified_user_genomes else 0

            len_mash_classified_ar53 = len(mash_classified_user_genomes['ar53']) \
                if 'ar53' in mash_classified_user_genomes else 0


            self.logger.info(f'{len_mash_classified_ar53+len_mash_classified_bac120} genome(s) have '
                             f'been classified using the ANI pre-screening step.')

        if skip_ani_screen and ani_summary_files is not None:
            # if the ani_Screen step was run, we need to load the results from the ani_summary_files
            # and use them to generate the final taxonomy file
            mash_classified_user_genomes = self.load_fastani_results_pre_pplacer(ani_summary_files)

        output_files = {}

        for marker_set_id in ('ar53', 'bac120'):
            warning_counter, prodigal_fail_counter = 0, 0
            if marker_set_id == 'ar53':
                marker_summary_fh = CopyNumberFileAR53(align_dir, prefix)
                marker_summary_fh.read()
                if os.path.isfile(os.path.join(align_dir,
                                               PATH_AR53_USER_MSA.format(prefix=prefix))):
                    user_msa_file = os.path.join(align_dir,
                                                 PATH_AR53_USER_MSA.format(prefix=prefix))
                else:
                    user_msa_file = os.path.join(align_dir,
                                                 PATH_AR53_USER_MSA.format(prefix=prefix) + '.gz')
                summary_file = ClassifySummaryFileAR53(out_dir, prefix)
                red_dict_file = REDDictFileAR53(out_dir, prefix)
                disappearing_genomes_file = DisappearingGenomesFileAR53(out_dir, prefix)
                pplacer_classify_file = PplacerClassifyFileAR53(out_dir, prefix)
            elif marker_set_id == 'bac120':
                marker_summary_fh = CopyNumberFileBAC120(align_dir, prefix)
                marker_summary_fh.read()
                if os.path.isfile(os.path.join(align_dir,
                                               PATH_BAC120_USER_MSA.format(prefix=prefix))):
                    user_msa_file = os.path.join(align_dir,
                                                 PATH_BAC120_USER_MSA.format(prefix=prefix))
                else:
                    user_msa_file = os.path.join(align_dir,
                                                 PATH_BAC120_USER_MSA.format(prefix=prefix) + '.gz')
                summary_file = ClassifySummaryFileBAC120(out_dir, prefix)
                red_dict_file = REDDictFileBAC120(out_dir, prefix)
                disappearing_genomes_file = DisappearingGenomesFileBAC120(out_dir, prefix)
                pplacer_classify_file = PplacerClassifyFileBAC120(out_dir, prefix)
                if not fulltreeopt:
                    tree_mapping_file = GenomeMappingFile(out_dir, prefix)
            else:
                raise GenomeMarkerSetUnknown('There was an error determining the marker set.')

            if (not os.path.exists(user_msa_file)) or (os.path.getsize(user_msa_file) < 30):
                # file will not exist if there are no User genomes from a given domain
                #
                # But if there is Unclassified genomes without domain,
                # they still have to be written in the bac120 summary file:
                if marker_set_id == 'bac120':
                    # Add failed genomes from prodigal and genomes with no markers in the bac120 summary file
                    # This is an executive direction: failed prodigal and genomes with no markers are not bacterial or archaeal
                    # but they need to be included in one of the summary file
                    prodigal_failed_counter = self.add_failed_genomes_to_summary(align_dir, summary_file, prefix)
                    if summary_file.has_row():
                        summary_file.write()
                        output_files.setdefault(marker_set_id, []).append(summary_file.path)
                        # Symlink to the summary file from the root
                        symlink_f(PATH_BAC120_SUMMARY_OUT.format(prefix=prefix),
                                  os.path.join(out_dir, os.path.basename(PATH_BAC120_SUMMARY_OUT.format(prefix=prefix))))
                        if prodigal_failed_counter > 0:
                            self.logger.warning(f"{prodigal_failed_counter} of {len(genomes)} "
                                                f"genome{'' if prodigal_failed_counter == 1 else 's'} "
                                                f"ha{'s' if prodigal_failed_counter == 1 else 've'} been labeled as 'Unclassified'.")

                continue

            msa_dict = read_fasta(user_msa_file)

            # Read the translation table summary file (identify).
            tln_table_summary_file = TlnTableSummaryFile(align_dir, prefix)
            tln_table_summary_file.read()

            percent_multihit_dict = self.parser_marker_summary_file(marker_summary_fh)

            # run pplacer to place bins in reference genome tree
            genomes_to_process = [seq_id for seq_id, _seq in read_seq(user_msa_file)]

            # if mash_classified_user_genomes is has key marker_set_id, we
            # need to add those genomes to the summary file and remove those genomes to the list of genomes
            # to process with pplacer
            if mash_classified_user_genomes and marker_set_id in mash_classified_user_genomes:
                list_summary_rows = mash_classified_user_genomes.get(marker_set_id)
                for row in list_summary_rows:
                    row,warning_counter = self._add_warning_to_row(row,
                                                   msa_dict,
                                                   tln_table_summary_file.genomes,
                                                   percent_multihit_dict,
                                                   warning_counter)
                    summary_file.add_row(row)
                    if row.gid in genomes_to_process:
                        genomes_to_process.remove(row.gid)

                # if genomes are classified with pre-screen, they need to be removed from the user_msa_file
                prescreened_msa_file_path = os.path.join(
                    out_dir, PATH_BAC120_PRESCREEN_MSA.format(prefix=prefix))
                #makes sure the path exists
                make_sure_path_exists(os.path.dirname(prescreened_msa_file_path))

                prescreened_msa_file = open(prescreened_msa_file_path, 'w')

                for gid in genomes_to_process:
                    prescreened_msa_file.write('>{}\n{}\n'.format(gid, msa_dict.get(gid)))
                prescreened_msa_file.close()
                user_msa_file = prescreened_msa_file_path




            # Write the RED dictionary to disk (intermediate file).
            red_dict_file.write()

            # if the new user_msa_file is empty, we need to skip the rest of the loop
            if not os.path.exists(user_msa_file) or os.path.getsize(user_msa_file) < 30:
                if summary_file.has_row():
                    summary_file.write()
                    # Symlink to the summary file from the root
                    if marker_set_id == 'bac120':
                        symlink_f(PATH_BAC120_SUMMARY_OUT.format(prefix=prefix),
                                  os.path.join(out_dir,
                                               os.path.basename(PATH_BAC120_SUMMARY_OUT.format(prefix=prefix))))
                    elif marker_set_id == 'ar53':
                        symlink_f(PATH_AR53_SUMMARY_OUT.format(prefix=prefix),
                                  os.path.join(out_dir, os.path.basename(PATH_AR53_SUMMARY_OUT.format(prefix=prefix))))
                continue


            if not fulltreeopt and marker_set_id == 'bac120':
                splitter = Split(self.order_rank, self.gtdb_taxonomy, self.reference_ids)

                debug_file = self._generate_summary_file(
                    marker_set_id, prefix, out_dir, debugopt, fulltreeopt)

                high_classify_tree = self.place_genomes(user_msa_file,
                                                        marker_set_id,
                                                        out_dir,
                                                        prefix,
                                                        scratch_dir,
                                                        'high')
                output_files.setdefault(marker_set_id, []).append(high_classify_tree)

                tree = self._assign_mrca_red(
                    high_classify_tree, marker_set_id, 'high')

                high_classification = splitter.get_high_pplacer_taxonomy(
                    out_dir, marker_set_id, prefix, user_msa_file, tree)

                disappearing_genomes = [seq_id for seq_id in genomes_to_process if seq_id not in high_classification]
                if disappearing_genomes:
                    for disappearing_genome in disappearing_genomes:
                        disappearing_genomes_file.add_genome(disappearing_genome, 'Backbone')

                tree_mapping_dict = {}
                tree_mapping_dict_reverse = {}
                with open(Config.CLASS_LEVEL_TREE_MAPPING_FILE) as ltmf:
                    for line in ltmf:
                        k, v = line.strip().split()
                        tree_mapping_dict[k] = v
                        tree_mapping_dict_reverse.setdefault(v, []).append(k)

                splitter = Split(self.order_rank, self.gtdb_taxonomy, self.reference_ids)
                sorted_high_taxonomy,warning_counter, len_sorted_genomes, high_taxonomy_used = splitter.map_high_taxonomy(
                    high_classification, tree_mapping_dict, summary_file, tree_mapping_file,msa_dict,
                    tln_table_summary_file.genomes,percent_multihit_dict,bac_ar_diff,warning_counter)

                if debugopt:
                    with open(out_dir + '/' + prefix + '_backbone_classification.txt', 'w') as htu:
                        for l, b in high_taxonomy_used.items():
                            htu.write(l + '\t' + str(b[0]) + '\t' + str(b[1]) + '\t' + str(b[2]) +'\n')

                self.logger.info(
                    f"{len_sorted_genomes} out of {len(genomes_to_process)} have an class assignments. Those genomes "
                    f"will be reclassified.")

                for idx, tree_iter in enumerate(
                        sorted(sorted_high_taxonomy, key=lambda z: len(sorted_high_taxonomy[z]), reverse=True)):
                    listg = sorted_high_taxonomy.get(tree_iter)
                    low_classify_tree, submsa_file_path = self._place_in_low_tree(
                        tree_iter, len(sorted_high_taxonomy), idx + 1, listg, msa_dict, marker_set_id, prefix,
                        scratch_dir, out_dir)
                    output_files.setdefault(marker_set_id, []).append(low_classify_tree)
                    genomes_to_process_subtree = [seq_id for seq_id, _seq in read_seq(submsa_file_path)]
                    mrca_lowtree = self._assign_mrca_red(
                        low_classify_tree, marker_set_id, 'low', tree_iter)
                    pplacer_classify_file = PplacerLowClassifyFileBAC120(out_dir, prefix, tree_iter)
                    pplacer_taxonomy_dict = self._get_pplacer_taxonomy(pplacer_classify_file,
                                                                       marker_set_id, user_msa_file, mrca_lowtree)
                    disappearing_genomes = [seq_id for seq_id in genomes_to_process_subtree if
                                            seq_id not in pplacer_taxonomy_dict]
                    if disappearing_genomes:
                        self.logger.warning(
                            f"{len(disappearing_genomes)} out of {len(genomes_to_process_subtree)} have not been"
                            f" properly placed by pplacer. This is a known issue with pplacer but we do "
                            f"not have a solution currently. Those missing genomes are written to "
                            f"the {disappearing_genomes_file.file_name} file. We recommend rerunning "
                            f"those genomes through GTDB-Tk to 'fix' the problem.")
                        for disappearing_genome in disappearing_genomes:
                            disappearing_genomes_file.add_genome(disappearing_genome, tree_iter)

                    class_level_classification, classified_user_genomes,warning_counter = self._parse_tree(mrca_lowtree, genomes, msa_dict,
                                                                  percent_multihit_dict, tln_table_summary_file.genomes,
                                                                  bac_ar_diff, submsa_file_path, red_dict_file.data,
                                                                  summary_file, pplacer_taxonomy_dict,warning_counter,
                                                                  high_classification, debug_file,skip_ani_screen, debugopt,
                                                                  tree_mapping_file, tree_iter,
                                                                  tree_mapping_dict_reverse)

                    if debugopt:
                        with open(out_dir + '/' + prefix + '_class_level_classification.txt', 'a') as olf:
                            for l, b in class_level_classification.items():
                                olf.write(l + '\t' + str(b) + '\n')
                            for l, b in classified_user_genomes.items():
                                olf.write(l + '\t' + str(b) + '\n')


                # add filtered genomes to the summary file
                warning_counter = self.add_filtered_genomes_to_summary(align_dir,warning_counter, summary_file, marker_set_id, prefix)

                # Add failed genomes from prodigal and genomes with no markers in the bac120 summary file
                # This is a executive direction: failed prodigal and genomes with
                # no markers are not bacterial or archaeal
                # but they need to be included in one of the summary file
                prodigal_failed_counter = self.add_failed_genomes_to_summary(align_dir, summary_file, prefix)
                warning_counter = warning_counter + prodigal_failed_counter

                # Symlink to the summary file from the root
                symlink_f(PATH_BAC120_SUMMARY_OUT.format(prefix=prefix),
                          os.path.join(out_dir, os.path.basename(PATH_BAC120_SUMMARY_OUT.format(prefix=prefix))))

                if debugopt:
                    debug_file.close()

            else:
                classify_tree = self.place_genomes(user_msa_file,
                                                   marker_set_id,
                                                   out_dir,
                                                   prefix,
                                                   scratch_dir)

                #genomes_to_process = [seq_id for seq_id, _seq in read_seq(user_msa_file)]

                # get taxonomic classification of each user genome
                debug_file = self._generate_summary_file(
                    marker_set_id, prefix, out_dir, debugopt, fulltreeopt)

                tree_to_process = self._assign_mrca_red(classify_tree, marker_set_id)
                pplacer_taxonomy_dict = self._get_pplacer_taxonomy(pplacer_classify_file,
                                                                   marker_set_id,
                                                                   user_msa_file,
                                                                   tree_to_process)

                disappearing_genomes = [seq_id for seq_id in genomes_to_process if seq_id not in pplacer_taxonomy_dict]
                class_level_classification, classified_user_genomes,warning_counter = self._parse_tree(tree_to_process, genomes, msa_dict, percent_multihit_dict,
                                 tln_table_summary_file.genomes,
                                 bac_ar_diff, user_msa_file, red_dict_file.data, summary_file,
                                 pplacer_taxonomy_dict,warning_counter, None,
                                 debug_file,skip_ani_screen, debugopt, None, None, None)
                # add filtered genomes to the summary file
                warning_counter = self.add_filtered_genomes_to_summary(align_dir,warning_counter, summary_file, marker_set_id, prefix)

                # Add failed genomes from prodigal and genomes with no markers in the bac120 summary file
                # This is a executive direction: failed prodigal and genomes with no markers are nit bacterial or archaeal
                # but they need to be included in one of the summary file
                if marker_set_id == 'bac120':
                    prodigal_fail_counter = self.add_failed_genomes_to_summary(align_dir, summary_file, prefix)
                    warning_counter = warning_counter + prodigal_fail_counter

                # Symlink to the summary file from the root
                if marker_set_id == 'bac120':
                    symlink_f(PATH_BAC120_SUMMARY_OUT.format(prefix=prefix),
                              os.path.join(out_dir, os.path.basename(PATH_BAC120_SUMMARY_OUT.format(prefix=prefix))))
                elif marker_set_id == 'ar53':
                    symlink_f(PATH_AR53_SUMMARY_OUT.format(prefix=prefix),
                              os.path.join(out_dir, os.path.basename(PATH_AR53_SUMMARY_OUT.format(prefix=prefix))))
                else:
                    raise GenomeMarkerSetUnknown('There was an error determining the marker set.')

                if disappearing_genomes:
                    for disappearing_genome in disappearing_genomes:
                        disappearing_genomes_file.add_genome(disappearing_genome, 'N/A')

                if debugopt:
                    debug_file.close()

            if not fulltreeopt and marker_set_id == 'bac120':
                tree_mapping_file.write()
                output_files.setdefault(marker_set_id, []).append(tree_mapping_file.path)


            if warning_counter > 0:
                sum_of_genomes =  len(genomes_to_process)+prodigal_failed_counter
                if mash_classified_user_genomes and marker_set_id in mash_classified_user_genomes:
                    sum_of_genomes +=len(mash_classified_user_genomes.get(marker_set_id))
                self.logger.warning(f"{warning_counter} of {sum_of_genomes} "
                                    f"genome{'' if warning_counter==1 else 's'}"
                                    f" ha{'s' if warning_counter==1 else 've'} a warning (see summary file).")

            # Write the summary file to disk.
            if disappearing_genomes_file.data:
                disappearing_genomes_file.write()
                output_files.setdefault(marker_set_id, []).append(disappearing_genomes_file.path)
            summary_file.write()
            output_files.setdefault(marker_set_id, []).append(summary_file.path)
        return output_files

    def _add_warning_to_row(self,row,msa_dict,
                            trans_table_dict,
                            percent_multihit_dict,
                            warning_counter):

        if row.gid in msa_dict:
            row.msa_percent = aa_percent_msa(msa_dict.get(row.gid))
            row.tln_table = trans_table_dict.get(row.gid)

        warnings = row.warnings
        if row.gid in percent_multihit_dict:
            if warnings is not None:
                warnings.append('Genome has more than {}% of markers with multiple hits'.format(
                    percent_multihit_dict.get(row.gid)))
            else:
                warnings = ['Genome has more than {}% of markers with multiple hits'.format(
                    percent_multihit_dict.get(row.gid))]
        if warnings:
            row.warnings = ';'.join(set(warnings))
            warning_counter += 1

        return row,warning_counter

    def _generate_summary_file(self, marker_set_id, prefix, out_dir, debugopt=None, fulltreeopt=None):
        debug_file = None

        if debugopt:
            debug_file = open(os.path.join(
                out_dir, prefix + '.{}.debug_file.tsv'.format(marker_set_id)), 'w')
            debug_file.write(
                "User genome\tRed value\tHigher rank\tHigher value\tLower rank\tLower value\tcase\tclosest_rank\ttool\n")
        return debug_file

    def _place_in_low_tree(self, tree_iter, number_low_trees, idx_tree, listg, msa_dict, marker_set_id, prefix,
                           scratch_dir, out_dir):
        make_sure_path_exists(os.path.join(
            out_dir, DIR_CLASS_LEVEL_PPLACER.format(iter=tree_iter)))
        submsa_file_path = os.path.join(
            out_dir, PATH_CLASS_LEVEL_BAC120_SUBMSA.format(iter=tree_iter))

        submsa_file = open(submsa_file_path, 'w')

        for gid in listg:
            submsa_file.write('>{}\n{}\n'.format(gid, msa_dict.get(gid)))
        submsa_file.close()
        low_classify_tree = self.place_genomes(submsa_file_path,
                                               marker_set_id,
                                               out_dir,
                                               prefix,
                                               scratch_dir,
                                               'low', tree_iter,
                                               number_low_trees, idx_tree)
        return low_classify_tree, submsa_file_path

    @staticmethod
    def _get_fastani_verification(tree, reference_ids, tt):
        """

        Parameters
        ----------
        tree : dendropy.Tree
            The input tree.

        reference_ids : frozenset
            A set containing the reference ID labels.

        Returns
        -------

        """
        """For each user genome, we select the first parent node with a label.
        If, while going up the tree, we find a node with only one reference
        genome, we select this reference genome as leaf_reference.
        """

        # Traverse up the tree, starting at each user leaf node.
        out = dict()
        number_comparison = 0
        qry_nodes = list(tree.leaf_node_iter(filter_fn=lambda x: x.taxon.label not in reference_ids))
        for leaf_node in tqdm_log(qry_nodes, unit='genome'):

            # Traverse up to find the first labelled parent node.
            par_node = leaf_node.parent_node
            leaf_ref_genome = None
            leaf_ref_genomes = [subnd for subnd in tt.get_leaf_nodes(par_node)
                                if subnd.taxon.label.replace("'", '') in reference_ids]
            if len(leaf_ref_genomes) == 1:
                leaf_ref_genome = leaf_ref_genomes[0]

            _support, parent_taxon, _aux_info = parse_label(par_node.label)
            # while par_node is not None and parent_taxon is empty,
            # we go up the tree
            while par_node is not None and not parent_taxon:
                par_node = par_node.parent_node
                if leaf_ref_genome is None:
                    leaf_ref_genomes = [subnd for subnd in tt.get_leaf_nodes(par_node)
                                        if subnd.taxon.label.replace("'", '') in reference_ids]
                    if len(leaf_ref_genomes) == 1:
                        leaf_ref_genome = leaf_ref_genomes[0]
                _support, parent_taxon, _aux_info = parse_label(par_node.label)

            # if the parent node is at the genus level
            parent_rank = parent_taxon.split(";")[-1]
            if parent_rank.startswith('g__'):
                # we get all the reference genomes under this genus
                list_subnode_initials = [subnd.taxon.label.replace("'", '')
                                         for subnd in tt.get_leaf_nodes(par_node)]
                if len(set(list_subnode_initials) & reference_ids) < 1:
                    raise GTDBTkExit(f"There are no reference genomes under '{parent_rank}'")
                else:
                    list_ref_genomes = [subnd for subnd in tt.get_leaf_nodes(par_node)
                                        if subnd.taxon.label.replace("'", '') in reference_ids]

                    # we pick the first 100 genomes closest (patristic distance) to the
                    # user genome under the same genus
                    dict_dist_refgenomes = calculate_patristic_distance(leaf_node,
                                                                        list_ref_genomes,
                                                                        tt)

                    sorted_l = sorted(iter(dict_dist_refgenomes.items()), key=itemgetter(1))
                    sorted_l = sorted_l[0:100]
                    number_comparison += len(sorted_l)
                    out[leaf_node] = {"potential_g": sorted_l,
                                      "pplacer_g": leaf_ref_genome}

            else:
                if leaf_ref_genome:
                    out[leaf_node] = {"potential_g": [(leaf_ref_genome, 0.0)],
                                      "pplacer_g": leaf_ref_genome}

        return out, qry_nodes

    def _classify_red_topology(self, tree, msa_dict, percent_multihit_dict, trans_table_dict, bac_ar_diff,
                               user_msa_file, red_dict,warning_counter, summary_file, pplacer_taxonomy_dict,
                               high_classification, debug_file, debugopt, classified_user_genomes,
                               unclassified_user_genomes, tt,tree_iter,tree_mapping_file, valid_classes,valid_phyla):
        user_genome_ids = set(read_fasta(user_msa_file).keys())
        user_genome_ids = user_genome_ids.difference(set(classified_user_genomes.keys()))
        class_level_classification = dict()
        for leaf in tree.leaf_node_iter(filter_fn=lambda x: x.taxon.label in user_genome_ids):

            # In some cases , pplacer can associate 2 user genomes
            # on the same parent node so we need to go up the tree
            # to find a node with a reference genome as leaf.
            cur_node = leaf.parent_node
            list_subnode = [subnd.taxon.label.replace("'", '')
                            for subnd in tt.get_leaf_nodes(cur_node)]
            while len(set(list_subnode) & self.reference_ids) < 1:
                cur_node = cur_node.parent_node
                list_subnode = [subnd.taxon.label.replace("'", '')
                                for subnd in tt.get_leaf_nodes(cur_node)]

            current_rel_list = cur_node.rel_dist

            parent_taxon_node = cur_node.parent_node
            _support, parent_taxon, _aux_info = parse_label(parent_taxon_node.label)

            while parent_taxon_node is not None and not parent_taxon:
                parent_taxon_node = parent_taxon_node.parent_node
                _support, parent_taxon, _aux_info = parse_label(parent_taxon_node.label)

            # is the node represent multiple ranks, we select the lowest one
            # i.e. if node is p__A;c__B;o__C we pick o__
            parent_rank = parent_taxon.split(";")[-1][0:3]
            parent_rel_dist = parent_taxon_node.rel_dist

            debug_info = [leaf.taxon.label, parent_rank, parent_rel_dist, '', '', '', '']

            child_taxons = []
            closest_rank = None
            detection = "taxonomic novelty determined using RED"
            # if the genome is not placed between the genus and
            # specie ranks
            if parent_rank != 'g__':
                # we select the child rank (if parent_rank = 'c__'
                # child rank will be 'o__)'
                child_rk = self.order_rank[self.order_rank.index(parent_rank) + 1]

                # get all reference genomes under the current node
                list_subnode = [childnd.taxon.label.replace("'", '')
                                for childnd in tt.get_leaf_nodes(cur_node)
                                if childnd.taxon.label in self.reference_ids]

                # get all names for the child rank
                list_ranks = [self.gtdb_taxonomy.get(name)[self.order_rank.index(child_rk)]
                              for name in list_subnode]

                # if there is just one rank name
                if len(set(list_ranks)) == 1:
                    for subranknd in cur_node.preorder_iter():
                        _support, subranknd_taxon, _aux_info = parse_label(subranknd.label)
                        if subranknd.is_internal() \
                                and subranknd_taxon is not None \
                                and subranknd_taxon.startswith(child_rk):
                            child_taxons = subranknd_taxon.split(";")
                            child_taxon_node = subranknd
                            child_rel_dist = child_taxon_node.rel_dist
                            break
                else:
                    # case 2a and 2b
                    closest_rank = parent_rank
                    detection = "taxonomic classification fully defined by topology"
            else:
                # case 1a
                closest_rank = parent_rank
                detection = "taxonomic classification fully defined by topology"

            # case 1b
            if len(child_taxons) == 0 and closest_rank is None:
                list_leaves = [childnd.taxon.label.replace("'", '')
                               for childnd in tt.get_leaf_nodes(cur_node)
                               if childnd.taxon.label in self.reference_ids]
                if len(list_leaves) != 1:
                    list_subrank = []
                    for leaf_subrank in list_leaves:
                        list_subrank.append(self.gtdb_taxonomy.get(leaf_subrank)
                                            [self.order_rank.index(parent_rank) + 1])
                    if len(set(list_subrank)) == 1:
                        print(leaf.taxon.label)
                        print(list_leaves)
                        print(list_subrank)
                        print(set(list_subrank))
                        raise GTDBTkExit('There should be only one leaf.')
                    else:
                        closest_rank = parent_rank
                        detection = "taxonomic classification fully defined by topology"
                list_leaf_ranks = self.gtdb_taxonomy.get(list_leaves[0])[
                                  self.order_rank.index(child_rk):-1]  # We remove the species name
                for leaf_taxon in reversed(list_leaf_ranks):
                    if leaf_taxon == list_leaf_ranks[0]:
                        if abs(current_rel_list - red_dict.get(leaf_taxon[:3])) < abs(
                                current_rel_list - red_dict.get(parent_rank)):
                            closest_rank = leaf_taxon[:3]
                            debug_info[3] = leaf_taxon
                            debug_info[5] = 'case 1b - III'
                            break
                    else:
                        pchildrank = list_leaf_ranks[list_leaf_ranks.index(leaf_taxon) - 1]
                        if abs(current_rel_list - red_dict.get(leaf_taxon[:3])) < abs(
                                current_rel_list - red_dict.get(pchildrank[:3])):
                            closest_rank = leaf_taxon[:3]
                            debug_info[1] = pchildrank
                            debug_info[2] = 1.0
                            debug_info[3] = leaf_taxon
                            debug_info[5] = 'case 1b - II'
                            break
                if closest_rank is None:
                    closest_rank = parent_rank
                    debug_info[3] = list_leaf_ranks[0]
                    debug_info[5] = 'case 1b - IV'

            # if there is multiple ranks on the child node (i.e genome between p__Nitrospirae and c__Nitrospiria;o__Nitrospirales;f__Nitropiraceae)
            # we loop through the list of rank from f_ to c_ rank
            for child_taxon in reversed(child_taxons):
                # if lower rank is c__Nitropiria
                if child_taxon == child_taxons[0]:
                    if (abs(current_rel_list - red_dict.get(child_taxon[:3])) < abs(
                            child_rel_dist - red_dict.get(child_taxon[:3])) and
                            abs(current_rel_list - red_dict.get(child_taxon[:3])) < abs(
                                current_rel_list - red_dict.get(parent_rank))):
                        debug_info[3] = ';'.join(child_taxons)
                        debug_info[4] = child_rel_dist
                        debug_info[5] = 'case 3b - II'
                        closest_rank = child_taxon[:3]
                    elif closest_rank is None:
                        closest_rank = parent_rank
                        debug_info[3] = ';'.join(child_taxons)
                        debug_info[4] = child_rel_dist
                        debug_info[5] = 'case 3b - III'
                else:
                    pchildrank = child_taxons[child_taxons.index(
                        child_taxon) - 1]
                    if (abs(current_rel_list - red_dict.get(child_taxon[:3])) < abs(
                            current_rel_list - red_dict.get(pchildrank[:3])) and
                            abs(current_rel_list - red_dict.get(child_taxon[:3])) < abs(
                                child_rel_dist - red_dict.get(child_taxon[:3]))):
                        closest_rank = child_taxon
                        debug_info[3] = ';'.join(child_taxons)
                        debug_info[4] = child_rel_dist
                        debug_info[5] = 'case 3b - I'
                        break

            # case 1b
            if closest_rank is None:
                raise Exception('closest rank is None')

            debug_info[6] = closest_rank

            list_subnode = [subnd.taxon.label.replace("'", '')
                            for subnd in tt.get_leaf_nodes(cur_node)]
            red_taxonomy = self._get_redtax(list_subnode, closest_rank)
            notes = []
            warnings = []

            class_level_tax = standardise_taxonomy(';'.join(red_taxonomy)).split(';')
            standardised_red_tax= standardise_taxonomy(';'.join(red_taxonomy))
            class_in_spe_tree = valid_classes

            class_level_classification[leaf.taxon.label] = standardised_red_tax

            if valid_classes:
                backbone_tax = high_classification.get(leaf.taxon.label).get('tk_tax_red').split(';')
                pplacer_taxonomy_to_report = ';'.join(class_level_tax)
                red_value_to_report = current_rel_list

                mapping_row = GenomeMappingFileRow()
                mapping_row.gid = leaf.taxon.label
                mapping_row.ani_classification = False
                mapping_row.mapped_tree = tree_iter



                if backbone_tax[self.CLASS_IDX] == 'c__' and class_level_tax[self.CLASS_IDX] == 'c__':
                    # no classification in class level tree so must
                    # use the classification in the backbone tree
                    final_split = backbone_tax
                    notes.append('classification based on placement in backbone tree')
                    mapping_row.mapped_tree = 'backbone'
                    mapping_row.rule = 'Rule 1'
                    pplacer_taxonomy_to_report = ';'.join(backbone_tax)
                    red_value_to_report = high_classification.get(leaf.taxon.label).get('rel_dist')

                elif class_level_tax[self.CLASS_IDX] in valid_classes \
                        and backbone_tax[self.CLASS_IDX] == class_level_tax[self.CLASS_IDX]:
                    # backbone and class-level trees have same class classification so we trust the
                    # class-level classification down to the rank of class
                    final_split = class_level_tax
                    mapping_row.rule = 'Rule 2'
                    notes.append('classification based on placement in class-level tree')

                elif class_level_tax[self.CLASS_IDX] in valid_classes:
                    # class classification in the class-level tree
                    # is to a valid class (i.e. a class fully contained
                    # in the class-level tree) so we trust the class-level
                    # classification despite the incongruency with the
                    # backbone tree
                    final_split = class_level_tax
                    mapping_row.rule = 'Rule 3'
                    notes.append('classification based on placement in class-level tree')
                    warnings.append('backbone tree has incongruent class classification')

                elif class_level_tax[self.PHYLUM_IDX] in valid_phyla \
                        and backbone_tax[self.PHYLUM_IDX] == class_level_tax[self.PHYLUM_IDX]:
                    # backbone and class-level trees have same
                    # phylum classification so we trust the
                    # class-level classification down to the
                    # rank of phylum
                    final_split = limit_rank(class_level_tax, self.PHYLUM_IDX)
                    mapping_row.rule = 'Rule 4'
                    notes.append('classification based on placement in class-level tree')
                    warnings.append('classification restricted to valid phylum')

                elif class_level_tax[self.PHYLUM_IDX] in valid_phyla:
                    # phylum classification in the class-level tree
                    # is to a valid phylum (i.e. a phylum from a valid class
                    # in the class-level tree) so we trust the class-level
                    # classification despite the incongruency with the
                    # backbone tree down to the rank of phylum
                    final_split = limit_rank(class_level_tax, self.PHYLUM_IDX)
                    mapping_row.rule = 'Rule 5'
                    notes.append('classification based on placement in class-level tree')
                    warnings.append('backbone tree has incongruent phylum classification')

                else:
                    # class-level classification is questionable since it is
                    # not to a valid taxa so use the consensus between
                    # the backbone and class-level classifications
                    consensus = []
                    for idx, (backbone_taxon, order_taxon) in enumerate(zip(backbone_tax, class_level_tax)):
                        if backbone_taxon == order_taxon:
                            consensus.append(backbone_taxon)
                        else:
                            consensus.append(self.order_rank[idx])
                    final_split = consensus
                    mapping_row.rule = 'Rule 6'
                    notes.append('classification based on consensus between backbone and class-level tree')

                summary_row = ClassifySummaryFileRow()
                if leaf.taxon.label in unclassified_user_genomes:
                    summary_row = unclassified_user_genomes.get(leaf.taxon.label)
                    if summary_row.note == '':
                        summary_row.note = None
                summary_row.gid = leaf.taxon.label
                summary_row.classification = ';'.join(final_split)
                summary_row.pplacer_tax = pplacer_taxonomy_to_report
                summary_row.red_value = red_value_to_report
                if summary_row.classification_method is None:
                    summary_row.classification_method = detection
                summary_row.msa_percent = aa_percent_msa(msa_dict.get(summary_row.gid))
                summary_row.tln_table = trans_table_dict.get(summary_row.gid)

                tree_mapping_file.add_row(mapping_row)

            else:
                del debug_info[0]
                summary_row = ClassifySummaryFileRow()
                if leaf.taxon.label in unclassified_user_genomes:
                    summary_row = unclassified_user_genomes.get(leaf.taxon.label)
                    if summary_row.note == '':
                        summary_row.note = None
                summary_row.gid = leaf.taxon.label
                summary_row.classification = standardise_taxonomy(';'.join(red_taxonomy))
                summary_row.pplacer_tax = pplacer_taxonomy_dict.get(leaf.taxon.label)
                if summary_row.classification_method is None:
                    summary_row.classification_method = detection
                summary_row.msa_percent = aa_percent_msa(msa_dict.get(summary_row.gid))
                summary_row.tln_table = trans_table_dict.get(summary_row.gid)
                summary_row.red_value = current_rel_list

            if summary_row.gid in percent_multihit_dict:
                warnings.append('Genome has more than {}% of markers with multiple hits'.format(
                    percent_multihit_dict.get(summary_row.gid)))
            if summary_row.gid in bac_ar_diff:
                warnings.append('Genome domain questionable ( {}% Bacterial, {}% Archaeal)'.format(
                    bac_ar_diff.get(summary_row.gid).get('bac120'),
                    bac_ar_diff.get(summary_row.gid).get('ar53')))

            if len(warnings) > 0:
                if summary_row.warnings is not None:
                    warnings.extend(summary_row.warnings.split(';'))
                summary_row.warnings = ';'.join(set(warnings))
                warning_counter += 1
            if len(notes) > 0:
                if summary_row.note is not None:
                    notes.extend(summary_row.note.split(';'))
                summary_row.note = ';'.join(set(notes))
            summary_file.add_row(summary_row)

            if debugopt:
                debug_file.write('{0}\t{1}\t{2}\t{3}\n'.format(
                    leaf.taxon.label, current_rel_list, '\t'.join(str(x) for x in debug_info), detection))


        return class_level_classification,warning_counter

    def _parse_tree(self, tree, genomes, msa_dict, percent_multihit_dict, trans_table_dict, bac_ar_diff,
                    user_msa_file, red_dict, summary_file, pplacer_taxonomy_dict,warning_counter, high_classification,
                    debug_file,prescreening, debugopt, tree_mapping_file, tree_iter, tree_mapping_dict_reverse):
        # Genomes can be classified by using FastANI or RED values
        # We go through all leaves of the tree. if the leaf is a user
        # genome we take its parent node and look at all the leaves
        # for this node.

        # Persist descendant information for efficient traversal.
        tt = TreeTraversal()

        self.logger.log(Config.LOG_TASK, 'Traversing tree to determine classification method.')
        fastani_verification, qury_nodes = self._get_fastani_verification(tree, self.reference_ids, tt)


        #DEBUG: Skip FastANI step
        #fastani_verification = {}

        # we run a fastani comparison for each user genomes against the
        # selected genomes in the same genus
        ##if not prescreening and len(fastani_verification) > 0:
        if len(fastani_verification) > 0:
            fastani = FastANI(cpus=self.cpus, force_single=True)
            d_ani_compare, d_paths = self._get_fastani_genome_path(
                fastani_verification, genomes)
            self.logger.log(Config.LOG_TASK,
                            f'Calculating average nucleotide identity using '
                            f'FastANI (v{fastani.version}).')
            all_fastani_dict = fastani.run(d_ani_compare, d_paths)
        else:
            all_fastani_dict = {}



        classified_user_genomes, unclassified_user_genomes,warning_counter = self._sort_fastani_results(
            fastani_verification, pplacer_taxonomy_dict, all_fastani_dict, msa_dict, percent_multihit_dict,
            trans_table_dict, bac_ar_diff,warning_counter, summary_file)
        #if not prescreening:
        self.logger.info(f'{len(classified_user_genomes):,} genome(s) have '
                             f'been classified using FastANI and pplacer.')

        if tree_mapping_file:
            for genome in classified_user_genomes.keys():
                mapping_row = GenomeMappingFileRow()
                mapping_row.gid = genome
                mapping_row.ani_classification = True
                mapping_row.mapped_tree = tree_iter
                mapping_row.rule = 'Rule 7'
                tree_mapping_file.add_row(mapping_row)

        # Iterate over each leaf node that was not classified with FastANI.
        class_in_spe_tree = None
        phyla_in_spe_tree = None
        if tree_mapping_dict_reverse:
            class_in_spe_tree = tree_mapping_dict_reverse.get(tree_iter)
            phyla_in_spe_tree = self.get_authorised_rank(class_in_spe_tree,self.PHYLUM_IDX)


        class_level_classification,warning_counter = self._classify_red_topology(tree, msa_dict, percent_multihit_dict,
                                                                 trans_table_dict, bac_ar_diff, user_msa_file,
                                                                 red_dict,warning_counter, summary_file,
                                                                 pplacer_taxonomy_dict, high_classification,
                                                                 debug_file, debugopt, classified_user_genomes,
                                                                 unclassified_user_genomes, tt,tree_iter,tree_mapping_file,
                                                                 class_in_spe_tree, phyla_in_spe_tree)

        return class_level_classification,classified_user_genomes,warning_counter

    def _assign_mrca_red(self, input_tree, marker_set_id, levelopt=None, tree_iter=None):
        """Parse the pplacer tree and write the partial taxonomy for each user genome based on their placements

        Parameters
        ----------
        input_tree : pplacer tree
        marker_set_id : bacterial or archeal id (bac120 or ar53)

        Returns
        -------
        tree: pplacer tree with RED value added to nodes of interest

        """

        self.logger.info('Calculating RED values based on reference tree.')
        tree = dendropy.Tree.get_from_path(input_tree,
                                           schema='newick',
                                           rooting='force-rooted',
                                           preserve_underscores=True)

        # Persist descendant information for efficient traversal.
        tt = TreeTraversal()

        if levelopt is None:
            red_file = Config.MRCA_RED_BAC120
        elif levelopt == 'high':
            red_file = Config.BACKBONE_RED_FILE
        elif levelopt == 'low':
            red_file = Config.CLASS_LEVEL_RED_FILE.format(iter=tree_iter)
        if marker_set_id == 'ar53':
            red_file = Config.MRCA_RED_AR53

        # create map from leave labels to tree nodes
        leaf_node_map = {}
        for leaf in tree.leaf_node_iter():
            leaf_node_map[leaf.taxon.label] = leaf

        # parse RED file and associate reference RED value to reference node in
        # the tree
        reference_nodes = set()
        with open(red_file) as rf:
            for line in rf:
                label_ids, red_value = line.strip().split('\t')
                labels = label_ids.split('|')
                if len(labels) == 2:
                    taxa = [leaf_node_map[label].taxon for label in labels]
                    node = tree.mrca(taxa=taxa)
                elif len(labels) == 1:
                    node = leaf_node_map[labels[0]]

                node.rel_dist = float(red_value)
                reference_nodes.add(node)

        # For all leaf nodes that are not reference genomes
        # We only give RED value to added nodes placed on a reference edge ( between a reference parent and a reference child)
        # The new red value for the pplacer node =
        # RED_parent + (RED_child -RED_parent) * ( (pplacer_disttoroot - parent_disttoroot) / (child_disttoroot - parent_disttoroot) )

        count = 0
        for nd in tree.leaf_nodes():
            if nd not in reference_nodes:
                nd.rel_dist = 1.0
                pplacer_node = nd
                pplacer_parent_node = pplacer_node.parent_node

                while not bool(tt.get_leaf_nodes(pplacer_node) & reference_nodes):
                    pplacer_node = pplacer_parent_node
                    pplacer_parent_node = pplacer_node.parent_node

                # perform level-order tree search to find first child
                # node that is part of the reference set
                for child in pplacer_node.levelorder_iter():
                    if child in reference_nodes:
                        child_node = child
                        break

                # find first parent node that is part of the reference set
                while not pplacer_parent_node in reference_nodes:
                    pplacer_parent_node = pplacer_parent_node.parent_node

                # we go up the tree until we reach pplacer_parent_node
                current_node = child_node.parent_node
                edge_length = child_node.edge_length
                on_pplacer_branch = False
                pplacer_edge_length = 0

                while current_node != pplacer_parent_node:
                    if on_pplacer_branch or current_node == pplacer_node:
                        on_pplacer_branch = True
                        pplacer_edge_length += current_node.edge_length
                    edge_length += current_node.edge_length
                    current_node = current_node.parent_node

                ratio = pplacer_edge_length / edge_length

                branch_rel_dist = child_node.rel_dist - pplacer_parent_node.rel_dist
                branch_rel_dist = pplacer_parent_node.rel_dist + branch_rel_dist * ratio

                pplacer_node.rel_dist = branch_rel_dist

        return tree

    def _get_pplacer_taxonomy(self, pplacer_classify_file, marker_set_id, user_msa_file, tree):
        """Parse the pplacer tree and write the partial taxonomy for each user genome based on their placements

        Parameters
        ----------
        pplacer_classify_file : output file object to write
        marker_set_id : bacterial or archaeal id (bac120 or ar53)
        user_msa_file : msa file listing all user genomes for a certain domain
        tree : pplacer tree including the user genomes

        Returns
        -------
        dictionary[genome_label]=pplacer_taxonomy

        """

        # We get the pplacer taxonomy for comparison
        user_genome_ids = set(read_fasta(user_msa_file).keys())
        for leaf in tree.leaf_node_iter():
            if leaf.taxon.label in user_genome_ids:
                taxa = []
                cur_node = leaf
                while cur_node.parent_node:
                    _support, taxon, _aux_info = parse_label(cur_node.label)
                    if taxon:
                        for t in taxon.split(';')[::-1]:
                            taxa.append(t.strip())
                    cur_node = cur_node.parent_node
                taxa_str = ';'.join(taxa[::-1])
                pplacer_classify_file.add_genome(leaf.taxon.label,
                                                 standardise_taxonomy(taxa_str, marker_set_id))
        pplacer_classify_file.write()

        return pplacer_classify_file.data

    def add_filtered_genomes_to_summary(self, align_dir,warning_counter, summary_file, marker_set_id,prefix):
        if marker_set_id == 'bac120':
            filtered_file = os.path.join(align_dir,PATH_BAC120_FILTERED_GENOMES.format(prefix=prefix))
            domain = 'Bacteria'
        else:
            filtered_file = os.path.join(align_dir,PATH_AR53_FILTERED_GENOMES.format(prefix=prefix))
            domain = 'Archaea'

        with open(filtered_file) as fin:
            for line in fin:
                infos = line.strip().split('\t')
                summary_row = ClassifySummaryFileRow()
                summary_row.gid = infos[0]
                summary_row.classification = f'Unclassified {domain}'
                summary_row.warnings = infos[1]
                summary_file.add_row(summary_row)
                warning_counter += 1
        return warning_counter

    def add_failed_genomes_to_summary(self, align_dir, summary_file, prefix):
        warning_counter = 0
        prodigal_failed_file = os.path.join(align_dir, PATH_FAILS.format(prefix=prefix))
        align_failed_file = os.path.join(align_dir, PATH_FAILED_ALIGN_GENOMES.format(prefix=prefix))
        for failfile in (prodigal_failed_file, align_failed_file):
            if os.path.exists(failfile):
                with open(failfile) as pff:
                    for line in pff:
                        infos = line.strip().split('\t')
                        summary_row = ClassifySummaryFileRow()
                        summary_row.gid = infos[0]
                        summary_row.classification = f'Unclassified'
                        summary_row.warnings = infos[1]
                        summary_file.add_row(summary_row)
                        warning_counter += 1
        return warning_counter

    @staticmethod
    def formatnote(sorted_dict,gtdb_taxonomy,species_radius, labels):
        """Format the note field by concatenating all information in a sorted dictionary

        Parameters
        ----------
        sorted_dict : sorted dictionary listing reference genomes, ani and alignment fraction for a specific user genome
                    (genomeid, {ani: value, af: value})
        labels : array of label that are removed from the note field

        Returns
        -------
        string
            note field

        """
        gtdb_taxonomy = {canonical_gid(k): v for k, v in gtdb_taxonomy.items()}
        note_list = []
        for element in sorted_dict:
            if element[0] not in labels:
                note_str = "{}, {}, {}, {}, {}".format(element[0],
                                                       gtdb_taxonomy.get(
                                                           add_ncbi_prefix(canonical_gid(element[0])))[6],
                                                       species_radius.get(
                                                           element[0]),
                                                       round(
                                                           element[1].get('ani'), 2),
                                                       round(element[1].get('af'),3))
                note_list.append(note_str)
        return note_list


    def _sort_fastani_results_pre_pplacer(self,fastani_results,bac_ar_diff):
        """ When run mash/FastANI on all genomes before using pplacer, we need to sort those results and store them for
        a later use

        Parameters
        ----------
        fastani_results : dictionary listing the fastani ANI for each user genomes against the potential genomes

        """
        classified_user_genomes = {}

        # sort the dictionary by ani then af
        for gid in fastani_results.keys():
            thresh_results = [(ref_gid, hit) for (ref_gid, hit) in fastani_results[gid].items() if
                              hit['af'] >= Config.AF_THRESHOLD and hit['ani'] >= self.gtdb_radii.get_rep_ani(
                                  canonical_gid(ref_gid))]
            closest = sorted(thresh_results, key=lambda x: (-x[1]['ani'], -x[1]['af']))
            if len(closest) > 0:
                summary_row = ClassifySummaryFileRow()
                summary_row.gid = gid
                summary_row.classification_method = 'ani_screen'
                if len(closest) > 1:
                    other_ref = '; '.join(self.formatnote(
                        closest,self.gtdb_taxonomy,self.species_radius, [gid]))
                    if len(other_ref) == 0:
                        summary_row.other_related_refs = None
                    else:
                        summary_row.other_related_refs = other_ref
                summary_row.note = 'classification based on ANI only'

                fastani_matching_reference = closest[0][0]
                current_ani = fastani_results.get(gid).get(
                    fastani_matching_reference).get('ani')
                current_af = fastani_results.get(gid).get(
                    fastani_matching_reference).get('af')

                summary_row.fastani_ref = fastani_matching_reference
                summary_row.fastani_ref_radius = str(
                    self.species_radius.get(fastani_matching_reference))
                summary_row.fastani_tax = ";".join(self.gtdb_taxonomy.get(
                    add_ncbi_prefix(fastani_matching_reference)))
                summary_row.fastani_ani = round(current_ani, 2)
                summary_row.fastani_af = round(current_af, 3)
                taxa_str = ";".join(self.gtdb_taxonomy.get(
                    add_ncbi_prefix(fastani_matching_reference)))
                summary_row.classification = standardise_taxonomy(
                    taxa_str)

                warnings = []
                if gid in bac_ar_diff:
                    warnings.append('Genome domain questionable ( {}% Bacterial, {}% Archaeal)'.format(
                        bac_ar_diff.get(gid).get('bac120'),
                        bac_ar_diff.get(gid).get('ar53')))

                if len(warnings) > 0:
                    summary_row.warnings = warnings

                if (taxa_str.split(';')[0]).split('__')[1] == 'Bacteria':
                    domain = 'bac120'
                else:
                    domain = 'ar53'
                classified_user_genomes.setdefault(domain, []).append(summary_row)

        return classified_user_genomes




    def _sort_fastani_results(self, fastani_verification, pplacer_taxonomy_dict,
                              all_fastani_dict, msa_dict, percent_multihit_dict,
                              trans_table_dict, bac_ar_diff,warning_counter, summary_file):
        """Format the note field by concatenating all information in a sorted dictionary

        Parameters
        ----------
        fastani_verification : dictionary listing the potential genomes associated with a user genome
        d[user_genome] = {"potential_g": [(potential_genome_in_same_genus,patristic distance)],
        "pplacer_g": genome_of_reference_selected_by_pplacer(if any)}
        all_fastani_dict : dictionary listing the fastani ANI for each user genomes against the potential genomes
        d[user_genome]={ref_genome1:{"af":af,"ani":ani},ref_genome2:{"af":af,"ani":ani}}
        summaryfout: output file

        Returns
        -------
        classified_user_genomes: list of genomes where FastANI and Placement in the reference tree have
        predicted a taxonomy
        unclassified_user_genomes: dictionary of genomes where FastANI and Placement in the reference tree have not
        predicted a taxonomy

        """
        classified_user_genomes = {}
        unclassified_user_genomes = {}

        for userleaf, potential_nodes in fastani_verification.items():

            summary_row = ClassifySummaryFileRow()

            warnings = []
            if userleaf.taxon.label in percent_multihit_dict:
                warnings.append('Genome has more than {}% of markers with multiple hits'.format(
                    percent_multihit_dict.get(userleaf.taxon.label)))
            if userleaf.taxon.label in bac_ar_diff:
                warnings.append('Genome domain questionable ( {}% Bacterial, {}% Archaeal)'.format(
                    bac_ar_diff.get(userleaf.taxon.label).get('bac120'),
                    bac_ar_diff.get(userleaf.taxon.label).get('ar53')))

            if potential_nodes.get("pplacer_g"):
                pplacer_leafnode = potential_nodes.get("pplacer_g").taxon.label
                if pplacer_leafnode[0:3] in ['RS_', 'GB_']:
                    pplacer_leafnode = pplacer_leafnode[3:]
                if userleaf.taxon.label in all_fastani_dict:
                    # import IPython; IPython.embed()
                    prefilter_af_reference_dictionary = {k: v for k, v in
                                                         all_fastani_dict.get(userleaf.taxon.label).items() if v.get(
                            'af') >= self.af_threshold}
                    sorted_prefilter_af_dict = sorted(iter(prefilter_af_reference_dictionary.items()),
                                                      key=lambda _x_y1: (_x_y1[1]['ani'], _x_y1[1]['af']), reverse=True)

                    sorted_dict = sorted(iter(all_fastani_dict.get(
                        userleaf.taxon.label).items()), key=lambda _x_y: (_x_y[1]['ani'], _x_y[1]['af']), reverse=True)

                    fastani_matching_reference = None
                    if len(sorted_prefilter_af_dict) > 0:
                        if sorted_prefilter_af_dict[0][1].get('ani') >= self.species_radius.get(
                                sorted_prefilter_af_dict[0][0]):
                            fastani_matching_reference = sorted_prefilter_af_dict[0][0]
                            current_ani = all_fastani_dict.get(userleaf.taxon.label).get(
                                fastani_matching_reference).get('ani')
                            current_af = all_fastani_dict.get(userleaf.taxon.label).get(
                                fastani_matching_reference).get('af')
                        else:
                            warnings.append(
                                "Genome not assigned to closest species as it falls outside its pre-defined ANI radius")

                    taxa_str = ";".join(self.gtdb_taxonomy.get(
                        add_ncbi_prefix(pplacer_leafnode)))

                    summary_row.gid = userleaf.taxon.label

                    summary_row.pplacer_tax = pplacer_taxonomy_dict.get(userleaf.taxon.label)
                    summary_row.classification_method = 'taxonomic classification defined by topology and ANI'
                    summary_row.msa_percent = aa_percent_msa(
                        msa_dict.get(summary_row.gid))
                    summary_row.tln_table = trans_table_dict.get(summary_row.gid)
                    if len(warnings) > 0:
                        summary_row.warnings = ';'.join(warnings)

                    if fastani_matching_reference is not None:
                        summary_row.fastani_ref = fastani_matching_reference
                        summary_row.fastani_ref_radius = str(
                            self.species_radius.get(fastani_matching_reference))
                        summary_row.fastani_tax = ";".join(self.gtdb_taxonomy.get(
                            add_ncbi_prefix(fastani_matching_reference)))
                        summary_row.fastani_ani = round(current_ani, 2)
                        summary_row.fastani_af = round(current_af,3)
                        if pplacer_leafnode == fastani_matching_reference:
                            if taxa_str.endswith("s__"):
                                taxa_str = taxa_str + pplacer_leafnode
                            summary_row.classification = standardise_taxonomy(
                                taxa_str)
                            summary_row.closest_placement_ref = summary_row.fastani_ref
                            summary_row.closest_placement_radius = summary_row.fastani_ref_radius
                            summary_row.closest_placement_tax = summary_row.fastani_tax
                            summary_row.closest_placement_ani = summary_row.fastani_ani
                            summary_row.closest_placement_af = summary_row.fastani_af
                            summary_row.note = 'topological placement and ANI have congruent species assignments'
                            if len(sorted_dict) > 0:
                                other_ref = '; '.join(self.formatnote(
                                    sorted_dict,self.gtdb_taxonomy,self.species_radius, [fastani_matching_reference]))
                                if len(other_ref) == 0:
                                    summary_row.other_related_refs = None
                                else:
                                    summary_row.other_related_refs = other_ref

                        else:
                            taxa_str = ";".join(self.gtdb_taxonomy.get(
                                add_ncbi_prefix(fastani_matching_reference)))
                            summary_row.classification = standardise_taxonomy(
                                taxa_str)
                            summary_row.closest_placement_ref = pplacer_leafnode
                            summary_row.closest_placement_radius = str(
                                self.species_radius.get(pplacer_leafnode))
                            summary_row.closest_placement_tax = ";".join(self.gtdb_taxonomy.get(
                                add_ncbi_prefix(pplacer_leafnode)))
                            if pplacer_leafnode in all_fastani_dict.get(userleaf.taxon.label):
                                summary_row.closest_placement_ani = round(all_fastani_dict.get(
                                    userleaf.taxon.label).get(pplacer_leafnode).get('ani'), 2)
                                summary_row.closest_placement_af = round(all_fastani_dict.get(
                                    userleaf.taxon.label).get(pplacer_leafnode).get('af'),3)
                            summary_row.classification_method = 'ANI'

                            if len(sorted_dict) > 0:
                                other_ref = '; '.join(self.formatnote(
                                    sorted_dict,self.gtdb_taxonomy,self.species_radius, [fastani_matching_reference, pplacer_leafnode]))
                                if len(other_ref) == 0:
                                    summary_row.other_related_refs = None
                                else:
                                    summary_row.other_related_refs = other_ref
                        summary_file.add_row(summary_row)
                        classified_user_genomes[userleaf.taxon.label] = standardise_taxonomy(taxa_str)
                    else:
                        summary_row.closest_placement_ref = pplacer_leafnode
                        summary_row.closest_placement_radius = str(
                            self.species_radius.get(pplacer_leafnode))
                        summary_row.closest_placement_tax = ";".join(self.gtdb_taxonomy.get(
                            add_ncbi_prefix(pplacer_leafnode)))
                        if pplacer_leafnode in all_fastani_dict.get(userleaf.taxon.label):
                            summary_row.closest_placement_ani = round(all_fastani_dict.get(
                                userleaf.taxon.label).get(pplacer_leafnode).get('ani'), 2)
                            summary_row.closest_placement_af = round(all_fastani_dict.get(
                                userleaf.taxon.label).get(pplacer_leafnode).get('af'),3)

                        if len(sorted_dict) > 0:
                            other_ref = '; '.join(self.formatnote(
                                sorted_dict,self.gtdb_taxonomy,self.species_radius, [pplacer_leafnode]))
                            if len(other_ref) == 0:
                                summary_row.other_related_refs = None
                            else:
                                summary_row.other_related_refs = other_ref
                        unclassified_user_genomes[userleaf.taxon.label] = summary_row

            elif userleaf.taxon.label in all_fastani_dict:
                prefilter_af_reference_dictionary = {k: v for k, v in
                                                     all_fastani_dict.get(userleaf.taxon.label).items() if v.get(
                        'af') >= self.af_threshold}
                sorted_prefilter_af_dict = sorted(iter(prefilter_af_reference_dictionary.items()),
                                                  key=lambda _x_y1: (_x_y1[1]['ani'], _x_y1[1]['af']), reverse=True)
                sorted_dict = sorted(iter(all_fastani_dict.get(
                    userleaf.taxon.label).items()), key=lambda _x_y2: (_x_y2[1]['ani'], _x_y2[1]['af']), reverse=True)

                summary_row.gid = userleaf.taxon.label
                summary_row.pplacer_tax = pplacer_taxonomy_dict.get(
                    userleaf.taxon.label)
                summary_row.classification_method = 'ANI'
                summary_row.msa_percent = aa_percent_msa(
                    msa_dict.get(summary_row.gid))
                summary_row.tln_table = trans_table_dict.get(summary_row.gid)

                exception_genomes = []
                if len(sorted_prefilter_af_dict) > 0:

                    if len(sorted_dict) > 0:
                        other_ref = '; '.join(self.formatnote(
                            sorted_dict,self.gtdb_taxonomy,self.species_radius, exception_genomes))
                        if len(other_ref) == 0:
                            summary_row.other_related_refs = None
                        else:
                            summary_row.other_related_refs = other_ref

                    if len(warnings) > 0:
                        summary_row.warnings = ';'.join(warnings)
                        warning_counter += 1
                    if sorted_prefilter_af_dict[0][1].get('ani') >= self.species_radius.get(
                            sorted_prefilter_af_dict[0][0]):
                        fastani_matching_reference = sorted_prefilter_af_dict[0][0]
                        exception_genomes.append(fastani_matching_reference)

                        taxa_str = ";".join(self.gtdb_taxonomy.get(
                            add_ncbi_prefix(fastani_matching_reference)))
                        summary_row.classification = standardise_taxonomy(
                            taxa_str)

                        summary_row.fastani_ref = fastani_matching_reference
                        summary_row.fastani_ref_radius = str(
                            self.species_radius.get(fastani_matching_reference))
                        summary_row.fastani_tax = ";".join(self.gtdb_taxonomy.get(
                            add_ncbi_prefix(fastani_matching_reference)))
                        current_ani = all_fastani_dict.get(userleaf.taxon.label).get(
                            fastani_matching_reference).get('ani')
                        summary_row.fastani_ani = round(current_ani, 2)
                        current_af = all_fastani_dict.get(userleaf.taxon.label).get(
                            fastani_matching_reference).get('af')
                        summary_row.fastani_af = round(current_af,3)
                        summary_row.note = 'topological placement and ANI have incongruent species assignments'
                        if len(warnings) > 0:
                            summary_row.warnings = ';'.join(warnings)

                        summary_file.add_row(summary_row)
                        classified_user_genomes[userleaf.taxon.label] = standardise_taxonomy(taxa_str)
                    else:
                        warnings.append("Genome not assigned to closest species as "
                                     "it falls outside its pre-defined ANI radius")
                        summary_row.warnings = ';'.join(warnings)
                        summary_row.classification_method = 'taxonomic classification defined by topology and ANI'
                        unclassified_user_genomes[userleaf.taxon.label] = summary_row

                else:
                    if len(sorted_dict) > 0:
                        other_ref = '; '.join(self.formatnote(
                            sorted_dict,self.gtdb_taxonomy,self.species_radius, exception_genomes))
                        if len(other_ref) == 0:
                            summary_row.other_related_refs = None
                        else:
                            summary_row.other_related_refs = other_ref
                    unclassified_user_genomes[userleaf.taxon.label] = summary_row
        return classified_user_genomes, unclassified_user_genomes,warning_counter

    def _get_redtax(self, list_subnode, closest_rank):
        """
        Provide a taxonomy string to a user genome based on the reference genomes of the same clade.
        If the clade contains multiple reference genomes we are comparing their taxonomies.
        -If all reference genomes have the same taxonomy up to the 'closest rank' ,
        the taxonomy string including the closest rank is returned.
        -If **NOT** all reference genomes have the same taxonomy up to the 'closest rank',
        the taxonomy string **NOT** including the closest rank is returned.

        Parameters
        ----------
        list_subnode : list of leaf nodes including multiple reference genome.
        closest_rank : last rank of the reference taxonomy

        Returns
        -------
        string
            Taxonomy string.

        """

        subtax, multirefrank = self._parse_subnodes(list_subnode, closest_rank)
        # if all orders in the list are the same, the user genomes gets the
        # same order
        if len(set(multirefrank)) == 1:
            # case d
            subtax.append(multirefrank[0])
        else:
            # otherwise it's stored as undefined
            # case a,b
            subtax.append(closest_rank + "undefined")
        return subtax

    def _parse_subnodes(self, list_subnode, closest_rank):
        subtax = []
        multirefrank = []
        initial_loop = True
        for item in list_subnode:
            # We get the taxonomy of all reference genomes
            if item in self.reference_ids:
                taxonomy_from_file = self.gtdb_taxonomy.get(item)
                # we store the selected rank (i.e. order) for each reference
                # genome
                for rank in taxonomy_from_file:
                    if rank.startswith(closest_rank):
                        multirefrank.append(rank)
                        initial_loop = False
                        break
                    elif initial_loop:
                        # The first iteration is used to stored upper level (
                        # i.e. domain,phylum,class )
                        subtax.append(rank)
        return subtax, multirefrank

    def _calculate_red_distances(self, input_tree, out_dir):
        """
        Provide a taxonomy string to a user genome based on the reference genomes of the same clade.
        If the clade contains multiple reference genomes we are comparing their taxonomies.
        -If all reference genomes have the same taxonomy up to the 'closest rank' ,
        the taxonomy string including the closest rank is returned.
        -If **NOT** all reference genomes have the same taxonomy up to the 'closest rank',
        the taxonomy string **NOT** including the closest rank is returned.

        Parameters
        ----------
        list_subnode : list of leaf nodes including multiple reference genome.
        closest_rank : last rank of the reference taxonomy

        Returns
        -------
        dendropy.Tree
            Taxonomy string.
        """

        # read tree
        self.logger.info('Reading tree.')
        tree = dendropy.Tree.get_from_path(input_tree,
                                           schema='newick',
                                           rooting='force-rooted',
                                           preserve_underscores=True)

        self.logger.info('Reading taxonomy from file.')
        taxonomy = Taxonomy().read(Config.TAXONOMY_FILE)

        # determine taxa to be used for inferring distribution
        trusted_taxa = None
        taxa_for_dist_inference = self._filter_taxa_for_dist_inference(tree,
                                                                       taxonomy,
                                                                       trusted_taxa,
                                                                       Config.RED_MIN_CHILDREN,
                                                                       Config.RED_MIN_SUPPORT)

        phylum_rel_dists, rel_node_dists = self.median_rd_over_phyla(tree,
                                                                     taxa_for_dist_inference,
                                                                     taxonomy)

        # set edge lengths to median value over all rootings
        tree.seed_node.rel_dist = 0.0
        for n in tree.preorder_node_iter(lambda x: x != tree.seed_node):
            n.rel_dist = np_median(rel_node_dists[n.id])
            rd_to_parent = n.rel_dist - n.parent_node.rel_dist
            if rd_to_parent < 0:
                # This can occur since we are setting all nodes
                # to their median RED value.
                # self.logger.warning('Not all branches are positive after scaling.')
                pass
            n.edge_length = rd_to_parent

        if False:
            # These plots can be useful for debugging and internal use,
            # but are likely to be confusing to users.
            rd = RelativeDistance()

            input_tree_name = os.path.splitext(os.path.basename(input_tree))[0]
            plot_file = os.path.join(out_dir, '{}.png'.format(input_tree_name))
            rd._distribution_summary_plot(
                phylum_rel_dists, taxa_for_dist_inference, plot_file)

            gtdb_parent_ranks = Taxonomy().parents(taxonomy)
            median_outlier_table = os.path.join(
                out_dir, '{}.tsv'.format(input_tree_name))
            median_rank_file = os.path.join(
                out_dir, '{}.dict'.format(input_tree_name))
            rd._median_summary_outlier_file(phylum_rel_dists,
                                            taxa_for_dist_inference,
                                            gtdb_parent_ranks,
                                            median_outlier_table,
                                            median_rank_file,
                                            False)

            input_tree_name = os.path.splitext(os.path.basename(input_tree))[0]
            output_tree = os.path.join(
                out_dir, '{}.scaled.tree'.format(input_tree_name))
            tree.write_to_path(output_tree,
                               schema='newick',
                               suppress_rooting=True,
                               unquoted_underscores=True)

        return tree

    def _filter_taxa_for_dist_inference(self, tree, taxonomy, trusted_taxa, min_children, min_support):
        """Determine taxa to use for inferring distribution of relative divergences.

        Parameters
        ----------
        tree : Dendropy Tree
            Phylogenetic tree.
        taxonomy : d[taxon ID] -> [d__x; p__y; ...]
            Taxonomy for each taxon.
        trusted_taxa : iterable
            Trusted taxa to consider when inferring distribution.
        min_children : int
            Only consider taxa with at least the specified number of children taxa when inferring distribution.
        min_support : float
            Only consider taxa with at least this level of support when inferring distribution.
        """

        # determine children taxa for each named group
        taxon_children = Taxonomy().taxon_children(taxonomy)

        # get all named groups
        taxa_for_dist_inference = set()
        for taxon_id, taxa in taxonomy.items():
            for taxon in taxa:
                taxa_for_dist_inference.add(taxon)

        # sanity check species names as these are a common problem
        species = set()
        for taxon_id, taxa in taxonomy.items():
            if len(taxa) > Taxonomy.rank_index['s__']:
                species_name = taxa[Taxonomy.rank_index['s__']]
                valid, error_msg = True, None
                if species_name != 's__':
                    valid, error_msg = Taxonomy().validate_species_name(
                        species_name, require_full=True, require_prefix=True)
                if not valid:
                    print('[Warning] Species name {} for {} is invalid: {}'.format(
                        species_name, taxon_id, error_msg))
                    continue

                species.add(species_name)

        # restrict taxa to those with a sufficient number of named children
        # Note: a taxonomic group with no children will not end up in the
        # taxon_children data structure so care must be taken when applying
        # this filtering criteria.
        if min_children > 0:
            valid_taxa = set()
            for taxon, children_taxa in taxon_children.items():
                if len(children_taxa) >= min_children:
                    valid_taxa.add(taxon)

            taxa_for_dist_inference.intersection_update(valid_taxa)

            # explicitly add in the species since they have no
            # children and thus be absent from the taxon_child dictionary
            taxa_for_dist_inference.update(species)

        # restrict taxa used for inferring distribution to those with
        # sufficient support
        if min_support > 0:
            for node in tree.preorder_node_iter():
                if not node.label or node.is_leaf():
                    continue

                # check for support value
                support, taxon_name, _auxiliary_info = parse_label(node.label)

                if not taxon_name:
                    continue

                if support and float(support) < min_support:
                    taxa_for_dist_inference.difference_update([taxon_name])
                elif not support and min_support > 0:
                    # no support value, so inform user if they were trying to
                    # filter on this property
                    print(
                        '[Error] Tree does not contain support values. As such, --min_support should be set to 0.')
                    continue

        # restrict taxa used for inferring distribution to the trusted set
        if trusted_taxa:
            taxa_for_dist_inference = trusted_taxa.intersection(
                taxa_for_dist_inference)

        return taxa_for_dist_inference

    def median_rd_over_phyla(self,
                             tree,
                             taxa_for_dist_inference,
                             taxonomy):
        """Calculate the median relative divergence over all phyla rootings.

        Parameters
        ----------
        tree : Tree
          Dendropy tree.
        taxa_for_dist_inference : set
          Taxa to use for inference relative divergence distributions.
        taxonomy : d[taxon_id] -> [d__, p__, ..., s__]
          Taxonomy of extant taxa.
        """

        # get list of phyla level lineages
        all_phyla = self._get_phyla_lineages(tree)
        self.logger.info('Identified %d phyla.' % len(all_phyla))

        phyla = [p for p in all_phyla if p in taxa_for_dist_inference]
        self.logger.info(
            'Using %d phyla as rootings for inferring RED distributions.' % len(phyla))
        if len(phyla) < 2:
            self.logger.error('Rescaling requires at least 2 valid phyla.')
            sys.exit(-1)

        # give each node a unique id
        for i, n in enumerate(tree.preorder_node_iter()):
            n.id = i

        # calculate relative divergence for tree rooted on each phylum
        phylum_rel_dists = {}
        rel_node_dists = defaultdict(list)
        rd = RelativeDistance()
        for p in phyla:
            phylum = p.replace('p__', '').replace(' ', '_').lower()
            status_msg = '==> Calculating information with rooting on {}.              '.format(
                phylum.capitalize())
            sys.stdout.write('\r{}'.format(status_msg))
            sys.stdout.flush()

            cur_tree = self.root_with_outgroup(tree, taxonomy, p)

            # calculate relative distance to taxa
            rel_dists = rd.rel_dist_to_named_clades(cur_tree)
            rel_dists.pop(0, None)  # remove results for Domain

            # remove named groups in outgroup
            children = Taxonomy().children(p, taxonomy)
            for r in rel_dists.keys():
                rel_dists[r].pop(p, None)

            for t in children:
                for r in rel_dists.keys():
                    rel_dists[r].pop(t, None)

            phylum_rel_dists[phylum] = rel_dists

            # calculate relative distance to all nodes
            rd.decorate_rel_dist(cur_tree)

            # determine which lineages represents the 'ingroup'
            ingroup_subtree = None
            for c in cur_tree.seed_node.child_node_iter():
                _support, taxon_name, _auxiliary_info = parse_label(c.label)
                if not taxon_name or p not in taxon_name:
                    ingroup_subtree = c
                    break

            # do a preorder traversal of 'ingroup' and record relative
            # divergence to nodes
            for n in ingroup_subtree.preorder_iter():
                rel_node_dists[n.id].append(n.rel_dist)

        sys.stdout.write(
            '==> Inference for RED distributions finished.                         ')
        sys.stdout.flush()
        sys.stdout.write('\n')

        return phylum_rel_dists, rel_node_dists

    def _get_phyla_lineages(self, tree):
        """Get list of phyla level lineages.

        Parameters
        ----------
        tree : Dendropy Tree
            Phylogenetic tree.

        Returns
        -------
        list
            List of phyla level lineages.
        """
        phyla = []
        for node in tree.preorder_node_iter():
            if not node.label or node.is_leaf():
                continue

            _support, taxon_name, _auxiliary_info = parse_label(node.label)
            if taxon_name:
                taxa = [x.strip() for x in taxon_name.split(';')]
                if taxa[-1].startswith('p__'):
                    phyla.append(taxa[-1])

        return phyla

    def root_with_outgroup(self, input_tree, taxonomy, outgroup_taxa):
        """Reroot the tree using the given outgroup.

        Parameters
        ----------
        input_tree : Dendropy Tree
          Tree to rerooted.
        taxonomy : dict
            Taxonomy for taxa.
        outgroup_taxa : iterable
          Labels of taxa in outgroup.

        Returns
        -------
        Dendropy Tree
            Deep-copy of original tree rerooted on outgroup.
        """

        new_tree = input_tree.clone()

        outgroup = set()
        for genome_id, taxa in taxonomy.items():
            if outgroup_taxa in taxa:
                outgroup.add(genome_id)

        outgroup_in_tree = set()
        ingroup_in_tree = set()
        for n in new_tree.leaf_node_iter():
            if n.taxon.label in outgroup:
                outgroup_in_tree.add(n.taxon)
            else:
                ingroup_in_tree.add(n)

        if len(outgroup_in_tree) == 0:
            self.logger.warning('No outgroup taxa identified in the tree.')
            self.logger.warning('Tree was not rerooted.')
            sys.exit(0)

        # There is a complication here. We wish to find the MRCA of the outgroup
        # taxa. Finding the MRCA requires a rooted tree and we have no guarantee
        # that the tree isn't currently rooted within the outgroup clade. There is
        # also no way to identify a node that is guaranteed to be outside the outgroup
        # clade. As such, the tree is randomly rooted on a leaf node not in the outgroup.
        # This random re-rooting is performed until the MRCA does not spans all taxa in
        # the tree.

        leaves_in_tree = sum([1 for _ in new_tree.leaf_node_iter()])
        while True:
            rnd_ingroup_leaf = random.sample(ingroup_in_tree, 1)[0]
            new_tree.reroot_at_edge(rnd_ingroup_leaf.edge,
                                    length1=0.5 * rnd_ingroup_leaf.edge_length,
                                    length2=0.5 * rnd_ingroup_leaf.edge_length)

            mrca = new_tree.mrca(taxa=outgroup_in_tree)
            leaves_in_mrca = sum([1 for _ in mrca.leaf_iter()])
            if leaves_in_mrca != leaves_in_tree:
                break

        if leaves_in_mrca == leaves_in_tree:
            self.logger.error('The MRCA spans all taxa in the tree.')
            self.logger.error(
                'This indicating the selected outgroup is likely polyphyletic in the current tree.')
            self.logger.error(
                'This should never occur. Please report this as a bug.')
            sys.exit(-1)

        if mrca.edge_length is None:
            # self.logger.info('Tree appears to already be rooted on this outgroup.')
            pass
        else:
            new_tree.reroot_at_edge(mrca.edge,
                                    length1=0.5 * mrca.edge_length,
                                    length2=0.5 * mrca.edge_length)

        return new_tree

    def _get_fastani_genome_path(self, fastani_verification, genomes):
        """Generates a queue of comparisons to be made and the paths to
        the corresponding genome id."""
        dict_compare, dict_paths = dict(), dict()

        for qry_node, qry_dict in fastani_verification.items():
            user_label = qry_node.taxon.label
            dict_paths[user_label] = genomes[user_label]
            dict_compare[user_label] = set()
            for node in qry_dict.get('potential_g'):
                leafnode = node[0]
                shortleaf = leafnode.taxon.label
                if leafnode.taxon.label.startswith('GB_') or leafnode.taxon.label.startswith('RS_'):
                    shortleaf = leafnode.taxon.label[3:]
                # TODEL UBA genomes
                if shortleaf.startswith("UBA"):
                    ref_path = os.path.join(
                        Config.FASTANI_GENOMES,
                        'UBA',
                        shortleaf + Config.FASTANI_GENOMES_EXT)
                else:
                    ref_path = os.path.join(
                        Config.FASTANI_GENOMES,
                        self.parse_leaf_to_dir_path(shortleaf),
                        shortleaf + Config.FASTANI_GENOMES_EXT)

                if not os.path.isfile(ref_path):
                    raise GTDBTkExit(f'Reference genome missing from FastANI database: {ref_path}')

                dict_compare[user_label].add(shortleaf)
                dict_paths[shortleaf] = ref_path

        return dict_compare, dict_paths

    def get_authorised_rank(self, order_in_spe_tree,rank_index):
        results = []
        for k, v in self.gtdb_taxonomy.items():
            if v[self.order_rank.index(self.rank_of_interest)] in order_in_spe_tree:
                results.append(v[rank_index])
        return list(set(results))

    def load_fastani_results_pre_pplacer(self,ani_summary_files):
        """Load the FastANI results for the genomes classified with ANI Screen."""
        fastani_results={}
        for domain,ani_summary_file_path in ani_summary_files.items():
            ani_summary_file = ANISummaryFile(ani_summary_file_path)
            fastani_results[domain]=ani_summary_file.read()

        classified_user_genomes = {}

        # sort the dictionary by ani then af
        for domain,fastani_results in fastani_results.items():
            for gid,ani_results in fastani_results.items():
                summary_row = ClassifySummaryFileRow()
                summary_row.gid = gid
                summary_row.classification_method = 'ani_screen'
                summary_row.note = 'classification based on ANI only'

                fastani_matching_reference = list(ani_results)[0]

                if 'other_refs' in fastani_results.get(gid).get(
                    fastani_matching_reference):
                    summary_row.other_related_refs = fastani_results.get(gid).get(
                    fastani_matching_reference).get('other_refs')
                current_ani = fastani_results.get(gid).get(
                        fastani_matching_reference).get('ani')
                current_af = fastani_results.get(gid).get(
                        fastani_matching_reference).get('af')

                summary_row.fastani_ref = fastani_matching_reference
                summary_row.fastani_ref_radius = str(
                        self.species_radius.get(fastani_matching_reference))
                summary_row.fastani_tax = fastani_results.get(gid).get(
                        fastani_matching_reference).get('taxonomy')
                summary_row.fastani_ani = round(current_ani, 2)
                summary_row.fastani_af = round(current_af,3)
                taxa_str = ";".join(self.gtdb_taxonomy.get(
                        add_ncbi_prefix(fastani_matching_reference)))
                summary_row.classification = standardise_taxonomy(
                        taxa_str)

                classified_user_genomes.setdefault(domain, []).append(summary_row)
        return classified_user_genomes

    def convert_rows_to_dict(self,list_of_summary_rows):
        """Convert a list of rows to a dictionary."""
        dict_of_rows = {}
        for domain,rows in list_of_summary_rows.items():
            for row in rows:
                if row.other_related_refs:
                    infos ={row.gid:{row.fastani_ref: {'ani': row.fastani_ani,
                                                 'af': row.fastani_af,
                                                 'other_refs': row.other_related_refs}}}
                else:
                    infos ={row.gid:{row.fastani_ref: {'ani': row.fastani_ani, 'af': row.fastani_af}}}
                dict_of_rows.setdefault(domain, {}).update(infos)

        return dict_of_rows



