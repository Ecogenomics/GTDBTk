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

import logging
import multiprocessing
import random
import shutil
import sys
import tempfile
from collections import defaultdict
from operator import itemgetter

import dendropy
from numpy import median as np_median

import config.config as Config
from biolib_lite.common import remove_extension, make_sure_path_exists
from biolib_lite.execute import check_dependencies
from biolib_lite.newick import parse_label
from biolib_lite.seq_io import read_seq, read_fasta
from biolib_lite.taxonomy import Taxonomy
from gtdbtk.config.output import *
from gtdbtk.exceptions import GenomeMarkerSetUnknown
from gtdbtk.external.pplacer import Pplacer
from gtdbtk.markers import Markers
from relative_distance import RelativeDistance
from tools import add_ncbi_prefix, splitchunks, symlink_f

sys.setrecursionlimit(15000)


class Classify(object):
    """Determine taxonomic classification of genomes by ML placement."""

    def __init__(self, cpus=1):
        """Initialize."""

        check_dependencies(['pplacer', 'guppy', 'fastANI'])

        self.taxonomy_file = Config.TAXONOMY_FILE
        self.af_threshold = Config.AF_THRESHOLD
        self.gtdb_taxonomy = Taxonomy().read(self.taxonomy_file)

        self.order_rank = ["d__", "p__", "c__", "o__", 'f__', 'g__', 's__']

        self.logger = logging.getLogger('timestamp')
        self.cpus = cpus

        self.species_radius = self.parse_radius_file()

    def parse_radius_file(self):
        results = {}
        with open(Config.RADII_FILE) as f:
            for line in f:
                infos = line.strip().split('\t')
                gid = infos[1]
                if infos[1].startswith('GB_') or infos[1].startswith('RS_'):
                    gid = gid[3:]
                results[gid] = float(infos[2])
        return results

    def place_genomes(self,
                      user_msa_file,
                      marker_set_id,
                      out_dir,
                      prefix,
                      scratch_dir=None):
        """Place genomes into reference tree using pplacer."""
        # rename user MSA file for compatibility with pplacer
        if not user_msa_file.endswith('.fasta'):
            if marker_set_id == 'bac120':
                t = PATH_BAC120_USER_MSA.format(prefix=prefix)
            elif marker_set_id == 'ar122':
                t = PATH_AR122_USER_MSA.format(prefix=prefix)
            else:
                self.logger.error('There was an error determining the marker set.')
                raise GenomeMarkerSetUnknown

            shutil.copyfile(user_msa_file, t)
            user_msa_file = t

        # run pplacer to place bins in reference genome tree
        num_genomes = sum([1 for _seq_id, _seq in read_seq(user_msa_file)])

        # check if a scratch file is to be created
        pplacer_mmap_file = None
        if scratch_dir:
            self.logger.info(
                'Using a scratch file for pplacer allocations. This decreases memory usage and performance.')
            pplacer_mmap_file = ' --mmap-file {}'.format(
                os.path.join(scratch_dir, prefix + ".pplacer.scratch"))

        # get path to pplacer reference package
        if marker_set_id == 'bac120':
            self.logger.info(
                'Placing {} bacterial genomes into reference tree with pplacer (be patient).'.format(num_genomes))
            pplacer_ref_pkg = os.path.join(
                Config.PPLACER_DIR, Config.PPLACER_BAC120_REF_PKG)
        elif marker_set_id == 'ar122':
            self.logger.info(
                'Placing {} archaeal genomes into reference tree with pplacer (be patient).'.format(num_genomes))
            pplacer_ref_pkg = os.path.join(
                Config.PPLACER_DIR, Config.PPLACER_AR122_REF_PKG)
        elif marker_set_id == 'rps23':
            self.logger.info(
                'Placing {} genomes into reference tree with pplacer (be patient).'.format(num_genomes))
            pplacer_ref_pkg = os.path.join(
                Config.PPLACER_DIR, Config.PPLACER_RPS23_REF_PKG)
        else:
            self.logger.error('Unknown marker set: {}'.format(marker_set_id))
            raise GenomeMarkerSetUnknown('Unknown marker set: {}'.format(marker_set_id))

        # create pplacer output directory
        pplacer_out_dir = os.path.join(out_dir, DIR_PPLACER)
        if not os.path.exists(pplacer_out_dir):
            os.makedirs(pplacer_out_dir)

        # run pplacer
        if marker_set_id == 'bac120':
            pplacer_out = os.path.join(out_dir, PATH_BAC120_PPLACER_OUT)
            pplacer_json_out = os.path.join(out_dir, PATH_BAC120_PPLACER_JSON)
        elif marker_set_id == 'ar122':
            pplacer_out = os.path.join(out_dir, PATH_AR122_PPLACER_OUT)
            pplacer_json_out = os.path.join(out_dir, PATH_AR122_PPLACER_JSON)
        else:
            self.logger.error('There was an error determining the marker set.')
            raise GenomeMarkerSetUnknown

        pplacer = Pplacer()
        pplacer.run(self.cpus, 'WAG', pplacer_ref_pkg, pplacer_json_out, user_msa_file, pplacer_out, pplacer_mmap_file)

        # extract tree
        if marker_set_id == 'bac120':
            tree_file = os.path.join(out_dir, PATH_BAC120_TREE_FILE.format(prefix=prefix))
        elif marker_set_id == 'ar122':
            tree_file = os.path.join(out_dir, PATH_AR122_TREE_FILE.format(prefix=prefix))
        else:
            self.logger.error('There was an error determining the marker set.')
            raise GenomeMarkerSetUnknown

        pplacer.tog(pplacer_json_out, tree_file)

        # Symlink to the tree summary file
        if marker_set_id == 'bac120':
            symlink_f(PATH_BAC120_TREE_FILE.format(prefix=prefix),
                      os.path.join(out_dir, os.path.basename(PATH_BAC120_TREE_FILE.format(prefix=prefix))))
        elif marker_set_id == 'ar122':
            symlink_f(PATH_AR122_TREE_FILE.format(prefix=prefix),
                      os.path.join(out_dir, os.path.basename(PATH_AR122_TREE_FILE.format(prefix=prefix))))
        else:
            self.logger.error('There was an error determining the marker set.')
            raise GenomeMarkerSetUnknown

        return tree_file

    def standardise_taxonomy(self, taxstring, marker_set=None):
        """Create a 7 rank taxonomy string from an incomplete taxonomy string

        Parameters
        ----------
        taxstring : str
            incomplete taxonomy string
        marker_set : str
            The marker set to use.

        Returns
        -------
        string
            7 rank taxonomy string.
        """
        # return taxstring
        taxlist = taxstring.split(";")
        if marker_set == 'bac120':
            taxlist.insert(0, 'd__Bacteria')
        if marker_set == 'ar122':
            taxlist.insert(0, 'd__Archaea')
        taxlist.extend(self.order_rank[len(taxlist):])
        new_taxstring = ";".join(taxlist)
        return new_taxstring

    def _write_red_dict(self, out_dir, prefix, marker_set_id):
        """Write the RED value for each rank to a file

        Parameters
        ----------
        out_dir : output directory
        prefix : desired prefix for output files
        marker_set_id : bacterial or archeal id (bac120 or ar122)

        Returns
        -------
        dictionary
            dictionary[rank_prefix] = red_value

        """

        marker_dict = {}
        if marker_set_id == 'bac120':
            marker_dict = Config.RED_DIST_BAC_DICT
            out_path = os.path.join(out_dir, PATH_BAC120_RED_DICT.format(prefix=prefix))
        elif marker_set_id == 'ar122':
            marker_dict = Config.RED_DIST_ARC_DICT
            out_path = os.path.join(out_dir, PATH_AR122_RED_DICT.format(prefix=prefix))
        else:
            self.logger.error('There was an error determining the marker set.')
            raise GenomeMarkerSetUnknown

        make_sure_path_exists(os.path.dirname(out_path))

        with open(out_path, 'w') as reddictfile:
            reddictfile.write('Phylum\t{}\n'.format(marker_dict.get('p__')))
            reddictfile.write('Class\t{}\n'.format(marker_dict.get('c__')))
            reddictfile.write('Order\t{}\n'.format(marker_dict.get('o__')))
            reddictfile.write('Family\t{}\n'.format(marker_dict.get('f__')))
            reddictfile.write('Genus\t{}\n'.format(marker_dict.get('g__')))

        return marker_dict

    def _parse_red_dict(self, red_dist_dict):
        results = {}
        for k, v in red_dist_dict.iteritems():
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

    def parser_marker_summary_file(self, marker_summary_file, marker_set_id):
        results = {}
        with open(marker_summary_file, 'r') as msf:
            msf.readline()
            for line in msf:
                infos = line.strip().split('\t')
                if marker_set_id == "bac120":
                    multi_hits_percent = (100 * float(infos[2])) / \
                                         Config.BAC_MARKER_COUNT
                elif marker_set_id == "ar122":
                    multi_hits_percent = (100 * float(infos[2])) / \
                                         Config.AR_MARKER_COUNT
                # print (marker_set_id, float(infos[3]), multi_hits_percent)
                if multi_hits_percent >= Config.DEFAULT_MULTIHIT_THRESHOLD:
                    results[infos[0]] = round(multi_hits_percent, 1)
        return results

    def parse_trans_table_file(self, trans_table_file):
        results = {}
        with open(trans_table_file, 'r') as msf:
            for line in msf:
                infos = line.strip().split('\t')
                results[infos[0]] = infos[1]
        return results

    def run(self,
            genomes,
            align_dir,
            out_dir,
            prefix,
            scratch_dir=None,
            recalculate_red=None,
            debugopt=False):
        """Classify genomes based on position in reference tree."""

        _bac_gids, _ar_gids, bac_ar_diff = Markers().genome_domain(align_dir, prefix)

        for marker_set_id in ('ar122', 'bac120'):

            if marker_set_id == 'ar122':
                marker_summary_file = os.path.join(align_dir, PATH_AR122_MARKER_SUMMARY.format(prefix=prefix))
                user_msa_file = os.path.join(align_dir, PATH_AR122_USER_MSA.format(prefix=prefix))
            elif marker_set_id == 'bac120':
                marker_summary_file = os.path.join(align_dir, PATH_BAC120_MARKER_SUMMARY.format(prefix=prefix))
                user_msa_file = os.path.join(align_dir, PATH_BAC120_USER_MSA.format(prefix=prefix))
            else:
                self.logger.error('There was an error determining the marker set.')
                raise GenomeMarkerSetUnknown

            if (not os.path.exists(user_msa_file)) or (os.path.getsize(user_msa_file) < 30):
                # file will not exist if there are no User genomes from a
                # given domain
                continue

            percent_multihit_dict = self.parser_marker_summary_file(
                marker_summary_file, marker_set_id)

            trans_table_file = os.path.join(align_dir, PATH_TLN_TABLE_SUMMARY.format(prefix=prefix))
            trans_table_dict = self.parse_trans_table_file(trans_table_file)

            msa_dict = read_fasta(user_msa_file)

            classify_tree = self.place_genomes(user_msa_file,
                                               marker_set_id,
                                               out_dir,
                                               prefix,
                                               scratch_dir)

            # get taxonomic classification of each user genome
            tree = dendropy.Tree.get_from_path(classify_tree,
                                               schema='newick',
                                               rooting='force-rooted',
                                               preserve_underscores=True)

            if marker_set_id == 'bac120':
                path_summary = os.path.join(out_dir, PATH_BAC120_SUMMARY_OUT.format(prefix=prefix))
            elif marker_set_id == 'ar122':
                path_summary = os.path.join(out_dir, PATH_AR122_SUMMARY_OUT.format(prefix=prefix))
            else:
                self.logger.error('There was an error determining the marker set.')
                raise GenomeMarkerSetUnknown

            pplacer_taxonomy_dict = self._get_pplacer_taxonomy(out_dir, prefix, marker_set_id, user_msa_file, tree)

            summaryfout = open(path_summary, 'w')
            if debugopt:
                debugfile = open(os.path.join(
                    out_dir, prefix + '.{}.debug_file.tsv'.format(marker_set_id)), 'w')

            marker_dict = self._write_red_dict(
                out_dir, prefix, marker_set_id)

            summaryfout.write(
                "user_genome\tclassification\tfastani_reference\tfastani_reference_radius\tfastani_taxonomy\tfastani_ani\tfastani_af\t" +
                "closest_placement_reference\tclosest_placement_taxonomy\tclosest_placement_ani\tclosest_placement_af\tpplacer_taxonomy\t" +
                "classification_method\tnote\tother_related_references(genome_id,species_name,radius,ANI,AF)\taa_percent\ttranslation_table\tred_value\twarnings\n")
            if debugopt:
                debugfile.write(
                    "User genome\tRed value\tHigher rank\tHigher value\tLower rank\tLower value\tcase\tclosest_rank\ttool\n")

            # Genomes can be classified by using FastANI or RED values
            # We go through all leaves of the tree. if the leaf is a user
            # genome we take its parent node and look at all the leaves
            # for this node.
            all_fastani_dict = {}
            self.logger.info(
                'Calculating average nucleotide identity using FastANI.')
            fastani_verification = {}
            number_comparison = 0
            for userleaf in tree.leaf_node_iter():
                # for each user genome, we select the first parent node with a label.
                # if, while going up the tree, we find a node with only one
                # reference genome, we select this reference genome as
                # leaf_reference.
                if userleaf.taxon.label[0:3] not in ['RS_', 'GB_', 'UBA']:

                    par_node = userleaf.parent_node
                    leaf_ref_genome = None
                    leaf_ref_genomes = [subnd for subnd in par_node.leaf_iter(
                    ) if subnd.taxon.label.replace("'", '')[0:3] in ['RS_', 'GB_', 'UBA']]
                    if len(leaf_ref_genomes) == 1:
                        leaf_ref_genome = leaf_ref_genomes[0]

                    _support, parent_taxon, _aux_info = parse_label(
                        par_node.label)
                    # while par_node is not None and parent_taxon is empty,
                    # we go up the tree
                    while par_node is not None and not parent_taxon:
                        par_node = par_node.parent_node
                        if leaf_ref_genome is None:
                            leaf_ref_genomes = [subnd for subnd in par_node.leaf_iter(
                            ) if subnd.taxon.label.replace("'", '')[0:3] in ['RS_', 'GB_', 'UBA']]
                            if len(leaf_ref_genomes) == 1:
                                leaf_ref_genome = leaf_ref_genomes[0]
                        _support, parent_taxon, _aux_info = parse_label(
                            par_node.label)

                    # if the parent node is at the genus level
                    parent_rank = parent_taxon.split(";")[-1]
                    if parent_rank.startswith('g__'):
                        # we get all the reference genomes under this genus
                        list_subnode_initials = [subnd.taxon.label.replace(
                            "'", '')[0:3] for subnd in par_node.leaf_iter()]
                        if (list_subnode_initials.count('RS_') + list_subnode_initials.count(
                                'GB_') + list_subnode_initials.count('UBA')) < 1:
                            raise Exception(
                                "There is no reference genomes under '{}'".format('parent_rank'))
                        else:
                            dict_dist_refgenomes = {}
                            list_ref_genomes = [subnd for subnd in par_node.leaf_iter(
                            ) if subnd.taxon.label.replace("'", '')[0:3] in ['RS_', 'GB_', 'UBA']]
                            # we pick the first 100 genomes closest (patristic distance) to the
                            # user genome under the same genus
                            for ref_genome in list_ref_genomes:
                                taxon_labels = [
                                    userleaf.taxon.label, ref_genome.taxon.label]
                                mrca = tree.mrca(taxon_labels=taxon_labels)
                                # the following command is faster than
                                # calculating the patristic distance
                                dict_dist_refgenomes[ref_genome] = (userleaf.distance_from_root(
                                ) - mrca.distance_from_root()) + (
                                                                               ref_genome.distance_from_root() - mrca.distance_from_root())
                            sorted_l = sorted(
                                dict_dist_refgenomes.iteritems(), key=itemgetter(1))
                            sorted_l = sorted_l[0:100]
                            number_comparison += len(sorted_l)
                            fastani_verification[userleaf] = {
                                "potential_g": sorted_l, "pplacer_g": leaf_ref_genome}
                    else:
                        if leaf_ref_genome:
                            fastani_verification[userleaf] = {"potential_g": [
                                (leaf_ref_genome, 0.0)], "pplacer_g": leaf_ref_genome}

            # we run a fastani comparison for each user genomes against the
            # selected genomes in the same genus
            manager = multiprocessing.Manager()
            out_q = manager.dict()
            procs = []
            nprocs = self.cpus

            if len(fastani_verification) > 0:
                for item in splitchunks(fastani_verification, nprocs):
                    p = multiprocessing.Process(
                        target=self._fastaniWorker,
                        args=(item, genomes, out_q))
                    procs.append(p)
                    p.start()

                # Wait for all worker processes to finish
                for p in procs:
                    p.join()
                    if p.exitcode == 1:
                        raise ValueError("Stop!!")

                all_fastani_dict = dict(out_q)

            classified_user_genomes, unclassified_user_genomes = self._sort_fastani_results(
                fastani_verification, pplacer_taxonomy_dict, all_fastani_dict, msa_dict, percent_multihit_dict,
                trans_table_dict, bac_ar_diff, summaryfout)

            self.logger.info('{0} genomes have been classify using FastANI and pplacer.'.format(
                len(classified_user_genomes)))

            # If Fastani can't select a taxonomy for a genome, we use RED
            # distances

            if recalculate_red:
                tree_to_process = self._calculate_red_distances(
                    classify_tree, out_dir)

            else:
                tree_to_process = self._assign_mrca_red(
                    classify_tree, marker_set_id)

            user_genome_ids = set(read_fasta(user_msa_file).keys())
            # we remove ids already classified with FastANI
            user_genome_ids = user_genome_ids.difference(
                set(classified_user_genomes))
            for leaf in tree_to_process.leaf_node_iter():
                if leaf.taxon.label in user_genome_ids:
                    # In some cases , pplacer can associate 2 user genomes
                    # on the same parent node so we need to go up the tree
                    # to find a node with a reference genome as leaf.
                    cur_node = leaf.parent_node
                    list_subnode_initials = [subnd.taxon.label.replace(
                        "'", '')[0:3] for subnd in cur_node.leaf_iter()]
                    while 'RS_' not in list_subnode_initials and 'GB_' not in list_subnode_initials and 'UBA' not in list_subnode_initials:
                        cur_node = cur_node.parent_node
                        list_subnode_initials = [subnd.taxon.label.replace(
                            "'", '')[0:3] for subnd in cur_node.leaf_iter()]

                    current_rel_list = cur_node.rel_dist

                    parent_taxon_node = cur_node.parent_node
                    _support, parent_taxon, _aux_info = parse_label(
                        parent_taxon_node.label)

                    while parent_taxon_node is not None and not parent_taxon:
                        parent_taxon_node = parent_taxon_node.parent_node
                        _support, parent_taxon, _aux_info = parse_label(
                            parent_taxon_node.label)

                    # is the node represent multiple ranks, we select the lowest one
                    # i.e. if node is p__A;c__B;o__C we pick o__
                    parent_rank = parent_taxon.split(";")[-1][0:3]
                    parent_rel_dist = parent_taxon_node.rel_dist

                    debug_info = [leaf.taxon.label, parent_rank,
                                  parent_rel_dist, '', '', '', '']

                    child_taxons = []
                    closest_rank = None
                    detection = "taxonomic novelty determined using RED"
                    # if the genome is not placed between the genus and
                    # specie ranks
                    if parent_rank != 'g__':
                        # we select the child rank (if parent_rank = 'c__'
                        # child rank will be 'o__)'
                        child_rk = self.order_rank[self.order_rank.index(
                            parent_rank) + 1]

                        # get all reference genomes under the current node
                        list_subnode = [childnd.taxon.label.replace("'", '') for childnd in cur_node.leaf_iter(
                        ) if childnd.taxon.label[0:3] in ['RS_', 'UBA', 'GB_']]

                        # get all names for the child rank
                        list_ranks = [self.gtdb_taxonomy.get(
                            name)[self.order_rank.index(child_rk)] for name in list_subnode]

                        # if there is just one rank name
                        if len(set(list_ranks)) == 1:
                            for subranknd in cur_node.preorder_iter():
                                _support, subranknd_taxon, _aux_info = parse_label(
                                    subranknd.label)
                                if subranknd.is_internal() and subranknd_taxon is not None and subranknd_taxon.startswith(
                                        child_rk):
                                    child_taxons = subranknd_taxon.split(
                                        ";")
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
                        list_leaves = [childnd.taxon.label.replace("'", '') for childnd in cur_node.leaf_iter(
                        ) if childnd.taxon.label[0:3] in ['RS_', 'UBA', 'GB_']]
                        if len(list_leaves) != 1:
                            list_subrank = []
                            for leaf_subrank in list_leaves:
                                list_subrank.append(self.gtdb_taxonomy.get(
                                    leaf_subrank)[self.order_rank.index(parent_rank) + 1])
                            if len(set(list_subrank)) == 1:
                                print(list_leaves)
                                print(list_subrank)
                                raise Exception(
                                    'There should be only one leaf.')
                                sys.exit(-1)
                            else:
                                closest_rank = parent_rank
                                detection = "taxonomic classification fully defined by topology"
                        list_leaf_ranks = self.gtdb_taxonomy.get(
                            list_leaves[0])[self.order_rank.index(child_rk):-1]  # We remove the species name
                        for leaf_taxon in reversed(list_leaf_ranks):
                            if leaf_taxon == list_leaf_ranks[0]:
                                if abs(current_rel_list - marker_dict.get(leaf_taxon[:3])) < abs(
                                        (current_rel_list) - marker_dict.get(parent_rank)):
                                    closest_rank = leaf_taxon[:3]
                                    debug_info[3] = leaf_taxon
                                    debug_info[5] = 'case 1b - III'
                                    break
                            else:
                                pchildrank = list_leaf_ranks[list_leaf_ranks.index(
                                    leaf_taxon) - 1]
                                if abs(current_rel_list - marker_dict.get(leaf_taxon[:3])) < abs(
                                        current_rel_list - marker_dict.get(pchildrank[:3])):
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
                            if (abs(current_rel_list - marker_dict.get(child_taxon[:3])) < abs(
                                    child_rel_dist - marker_dict.get(child_taxon[:3])) and
                                    abs(current_rel_list - marker_dict.get(child_taxon[:3])) < abs(
                                        current_rel_list - marker_dict.get(parent_rank))):
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
                            if (abs(current_rel_list - marker_dict.get(child_taxon[:3])) < abs(
                                    current_rel_list - marker_dict.get(pchildrank[:3])) and
                                    abs(current_rel_list - marker_dict.get(child_taxon[:3])) < abs(
                                        child_rel_dist - marker_dict.get(child_taxon[:3]))):
                                closest_rank = child_taxon
                                debug_info[3] = ';'.join(child_taxons)
                                debug_info[4] = child_rel_dist
                                debug_info[5] = 'case 3b - I'
                                break

                    # case 1b
                    if closest_rank is None:
                        raise Exception('closest rank is None')

                    debug_info[6] = closest_rank

                    list_subnode = [subnd.taxon.label.replace(
                        "'", '') for subnd in cur_node.leaf_iter()]
                    red_taxonomy = self._get_redtax(
                        list_subnode, closest_rank)

                    del debug_info[0]

                    summary_list = [None] * 19
                    if leaf.taxon.label in unclassified_user_genomes:
                        summary_list = unclassified_user_genomes.get(
                            leaf.taxon.label)
                        if summary_list[13] == '':
                            summary_list[13] = None
                    summary_list[0] = leaf.taxon.label
                    summary_list[1] = self.standardise_taxonomy(
                        red_taxonomy)
                    summary_list[11] = pplacer_taxonomy_dict.get(leaf.taxon.label)
                    summary_list[12] = 'Placement'
                    summary_list[13] = detection
                    summary_list[15] = self.aa_percent_msa(
                        msa_dict.get(summary_list[0]))
                    summary_list[16] = trans_table_dict.get(
                        summary_list[0])
                    summary_list[17] = current_rel_list

                    notes = []
                    if summary_list[0] in percent_multihit_dict:
                        notes.append('Genome has more than {}% of markers with multiple hits'.format(
                            percent_multihit_dict.get(summary_list[0])))
                    if summary_list[0] in bac_ar_diff:
                        notes.append('Genome domain questionable ( {}% Bacterial, {}% Archaeal)'.format(
                            bac_ar_diff.get(summary_list[0]).get('bac120'),
                            bac_ar_diff.get(summary_list[0]).get('ar122')))

                    if len(notes) > 0:
                        summary_list[18] = ';'.join(notes)
                    summaryfout.write("{0}\n".format(
                        '\t'.join(['N/A' if x is None else str(x) for x in summary_list])))
                    if debugopt:
                        debugfile.write('{0}\t{1}\t{2}\t{3}\n'.format(
                            leaf.taxon.label, current_rel_list, '\t'.join(str(x) for x in debug_info), detection))
                    else:
                        'debug false'

            summaryfout.close()

            # Symlink to the summary file from the root
            if marker_set_id == 'bac120':
                symlink_f(PATH_BAC120_SUMMARY_OUT.format(prefix=prefix),
                          os.path.join(out_dir, os.path.basename(PATH_BAC120_SUMMARY_OUT.format(prefix=prefix))))
            elif marker_set_id == 'ar122':
                symlink_f(PATH_AR122_SUMMARY_OUT.format(prefix=prefix),
                          os.path.join(out_dir, os.path.basename(PATH_AR122_SUMMARY_OUT.format(prefix=prefix))))
            else:
                self.logger.error('There was an error determining the marker set.')
                raise GenomeMarkerSetUnknown

            if debugopt:
                debugfile.close()

    def _assign_mrca_red(self, input_tree, marker_set_id):
        """Parse the pplacer tree and write the partial taxonomy for each user genome based on their placements

        Parameters
        ----------
        input_tree : pplacer tree
        marker_set_id : bacterial or archeal id (bac120 or ar122)

        Returns
        -------
        tree: pplacer tree with RED value added to nodes of interest

        """

        self.logger.info('Calculating RED values based on reference tree.')
        tree = dendropy.Tree.get_from_path(input_tree,
                                           schema='newick',
                                           rooting='force-rooted',
                                           preserve_underscores=True)

        red_file = Config.MRCA_RED_BAC120
        if marker_set_id == 'ar122':
            red_file = Config.MRCA_RED_AR122

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
        reference_pplacer_node = {}
        for nd in tree.leaf_nodes():
            if nd not in reference_nodes:
                nd.rel_dist = 1.0
                pplacer_node = nd
                pplacer_parent_node = pplacer_node.parent_node

                while not bool(set(pplacer_node.leaf_nodes()) & reference_nodes):
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

    def _get_pplacer_taxonomy(self, out_dir, prefix, marker_set_id, user_msa_file, tree):
        """Parse the pplacer tree and write the partial taxonomy for each user genome based on their placements

        Parameters
        ----------
        out_dir : output directory
        prefix : desired prefix for output files
        marker_set_id : bacterial or archeal id (bac120 or ar122)
        user_msa_file : msa file listing all user genomes for a certain domain
        tree : pplacer tree including the user genomes

        Returns
        -------
        dictionary[genome_label]=pplacer_taxonomy

        """

        out_root = os.path.join(out_dir, 'classify', 'intermediate_results')
        make_sure_path_exists(out_root)
        result = {}

        if marker_set_id == 'bac120':
            out_pplacer = os.path.join(out_dir, PATH_BAC120_PPLACER_CLASS.format(prefix=prefix))
        elif marker_set_id == 'ar122':
            out_pplacer = os.path.join(out_dir, PATH_AR122_PPLACER_CLASS.format(prefix=prefix))
        else:
            self.logger.error('There was an error determining the marker set.')
            raise GenomeMarkerSetUnknown

        # We get the pplacer taxonomy for comparison
        with open(out_pplacer, 'w') as pplaceout:
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
                    pplaceout.write('{}\t{}\n'.format(
                        leaf.taxon.label, self.standardise_taxonomy(taxa_str, marker_set_id)))
                    result[leaf.taxon.label] = self.standardise_taxonomy(taxa_str, marker_set_id)
        return result

    def _formatnote(self, sorted_dict, labels):
        """Format the note field by concatenating all information in a sorted dictionary

        Parameters
        ----------
        sorted_dict : sorted dictionary listing reference genomes, ani and alignement fraction for a specific user genome
                    (genomeid, {ani: value, af: value})
        labels : array of label that are removed from the note field

        Returns
        -------
        string
            note field

        """
        note_list = []
        for element in sorted_dict:
            if element[0] not in labels:
                note_str = "{}, {}, {}, {}, {}".format(element[0],
                                                       self.gtdb_taxonomy.get(add_ncbi_prefix(element[0]))[6],
                                                       self.species_radius.get(element[0]),
                                                       round(element[1].get('ani'), 2),
                                                       element[1].get('af'))
                note_list.append(note_str)
        return note_list

    def aa_percent_msa(self, aa_string):
        aa_len = sum([1 for c in aa_string if c.isalpha()])
        aa_perc = float(aa_len) / len(aa_string)
        return round(aa_perc * 100, 2)

    def _sort_fastani_results(self, fastani_verification, pplacer_taxonomy_dict,
                              all_fastani_dict, msa_dict, percent_multihit_dict,
                              trans_table_dict, bac_ar_diff, summaryfout):
        """Format the note field by concatenating all information in a sorted dictionary

        Parameters
        ----------
        fastani_verification : dictionary listing the potential genomes associated with a user genome d[user_genome] = {"potential_g": [
                                    (potential_genome_in_same_genus,patristic distance)], "pplacer_g": genome_of_reference_selected_by_pplacer(if any)}
        all_fastani_dict : dictionary listing the fastani ANI for each user genomes against the potential genomes d[user_genome]={ref_genome1:{"af":af,"ani":ani},ref_genome2:{"af":af,"ani":ani}}
        summaryfout: output file 

        Returns
        -------
        classified_user_genomes: list of genomes where FastANI and Placement in the reference tree have predicted a taxonomy
        unclassified_user_genomes: dictionary of genomes where FastANI and Placement in the reference tree have not  predicted a taxonomy

        """
        classified_user_genomes = []
        unclassified_user_genomes = {}
        for userleaf, potential_nodes in fastani_verification.iteritems():
            summary_list = [None] * 19

            notes = []
            if userleaf.taxon.label in percent_multihit_dict:
                notes.append('Genome has more than {}% of markers with multiple hits'.format(
                    percent_multihit_dict.get(userleaf.taxon.label)))
            if userleaf.taxon.label in bac_ar_diff:
                notes.append('Genome domain questionable ( {}% Bacterial, {}% Archaeal)'.format(
                    bac_ar_diff.get(userleaf.taxon.label).get('bac120'),
                    bac_ar_diff.get(userleaf.taxon.label).get('ar122')))
            if len(notes) > 0:
                summary_list[18] = ';'.join(notes)

            if potential_nodes.get("pplacer_g"):
                pplacer_leafnode = potential_nodes.get("pplacer_g").taxon.label
                if pplacer_leafnode[0:3] in ['RS_', 'GB_']:
                    pplacer_leafnode = pplacer_leafnode[3:]
                if userleaf.taxon.label in all_fastani_dict:
                    # import IPython; IPython.embed()
                    prefilter_reference_dictionary = {k: v for k, v in
                                                      all_fastani_dict.get(userleaf.taxon.label).iteritems() if (
                                                                  v.get('ani') >= self.species_radius.get(k) and v.get(
                                                              'af') >= self.af_threshold)}
                    sorted_dict = sorted(all_fastani_dict.get(
                        userleaf.taxon.label).iteritems(), key=lambda (_x, y): (y['ani'], y['af']), reverse=True)
                    sorted_prefilter_dict = sorted(prefilter_reference_dictionary.iteritems(),
                                                   key=lambda (_x, y): (y['ani'], y['af']), reverse=True)

                    fastani_matching_reference = None
                    if len(sorted_prefilter_dict) > 0:
                        fastani_matching_reference = sorted_prefilter_dict[0][0]
                        current_ani = all_fastani_dict.get(userleaf.taxon.label).get(
                            fastani_matching_reference).get('ani')
                        current_af = all_fastani_dict.get(userleaf.taxon.label).get(
                            fastani_matching_reference).get('af')

                    taxa_str = ";".join(self.gtdb_taxonomy.get(
                        add_ncbi_prefix(pplacer_leafnode)))

                    summary_list[0] = userleaf.taxon.label

                    summary_list[11] = pplacer_taxonomy_dict.get(userleaf.taxon.label)
                    summary_list[12] = 'ANI/Placement'
                    summary_list[15] = self.aa_percent_msa(
                        msa_dict.get(summary_list[0]))
                    summary_list[16] = trans_table_dict.get(summary_list[0])

                    if fastani_matching_reference is not None:
                        summary_list[2] = fastani_matching_reference
                        summary_list[3] = str(
                            self.species_radius.get(fastani_matching_reference))
                        summary_list[4] = ";".join(self.gtdb_taxonomy.get(
                            add_ncbi_prefix(fastani_matching_reference)))
                        summary_list[5] = round(current_ani, 2)
                        summary_list[6] = current_af
                        if pplacer_leafnode == fastani_matching_reference:
                            if taxa_str.endswith("s__"):
                                taxa_str = taxa_str + pplacer_leafnode
                            summary_list[1] = self.standardise_taxonomy(
                                taxa_str)
                            summary_list[7] = summary_list[2]
                            summary_list[8] = summary_list[4]
                            summary_list[9] = summary_list[5]
                            summary_list[10] = summary_list[6]
                            summary_list[13] = 'topological placement and ANI have congruent species assignments'
                            if len(sorted_dict) > 0:
                                other_ref = '; '.join(self._formatnote(
                                    sorted_dict, [fastani_matching_reference]))
                                if len(other_ref) == 0:
                                    summary_list[14] = None
                                else:
                                    summary_list[14] = other_ref

                        else:
                            taxa_str = ";".join(self.gtdb_taxonomy.get(
                                add_ncbi_prefix(fastani_matching_reference)))
                            summary_list[1] = self.standardise_taxonomy(
                                taxa_str)
                            summary_list[7] = pplacer_leafnode
                            summary_list[8] = ";".join(self.gtdb_taxonomy.get(
                                add_ncbi_prefix(pplacer_leafnode)))
                            if pplacer_leafnode in all_fastani_dict.get(userleaf.taxon.label):
                                summary_list[9] = round(all_fastani_dict.get(
                                    userleaf.taxon.label).get(pplacer_leafnode).get('ani'), 2)
                                summary_list[10] = all_fastani_dict.get(
                                    userleaf.taxon.label).get(pplacer_leafnode).get('af')
                            summary_list[13] = 'topological placement and ANI have incongruent species assignments'
                            summary_list[12] = 'ANI'

                            if len(sorted_dict) > 0:
                                other_ref = '; '.join(self._formatnote(
                                    sorted_dict, [fastani_matching_reference, pplacer_leafnode]))
                                if len(other_ref) == 0:
                                    summary_list[14] = None
                                else:
                                    summary_list[14] = other_ref

                        summaryfout.write("{}\n".format(
                            '\t'.join(['N/A' if x is None else str(x) for x in summary_list])))
                        classified_user_genomes.append(userleaf.taxon.label)
                    else:
                        summary_list[7] = pplacer_leafnode
                        summary_list[8] = ";".join(self.gtdb_taxonomy.get(
                            add_ncbi_prefix(pplacer_leafnode)))
                        if pplacer_leafnode in all_fastani_dict.get(userleaf.taxon.label):
                            summary_list[9] = round(all_fastani_dict.get(
                                userleaf.taxon.label).get(pplacer_leafnode).get('ani'), 2)
                            summary_list[10] = all_fastani_dict.get(
                                userleaf.taxon.label).get(pplacer_leafnode).get('af')

                        if len(sorted_dict) > 0:
                            other_ref = '; '.join(self._formatnote(
                                sorted_dict, [pplacer_leafnode]))
                            if len(other_ref) == 0:
                                summary_list[14] = None
                            else:
                                summary_list[14] = other_ref
                        unclassified_user_genomes[userleaf.taxon.label] = summary_list

            elif userleaf.taxon.label in all_fastani_dict:
                # import IPython; IPython.embed()
                prefilter_reference_dictionary = {k: v for k, v in
                                                  all_fastani_dict.get(userleaf.taxon.label).iteritems() if (
                                                              v.get('ani') >= self.species_radius.get(k) and v.get(
                                                          'af') >= self.af_threshold)}
                sorted_dict = sorted(all_fastani_dict.get(
                    userleaf.taxon.label).iteritems(), key=lambda (_x, y): (y['ani'], y['af']), reverse=True)
                sorted_prefilter_dict = sorted(prefilter_reference_dictionary.iteritems(),
                                               key=lambda (_x, y): (y['ani'], y['af']), reverse=True)

                summary_list[0] = userleaf.taxon.label
                summary_list[11] = pplacer_taxonomy_dict.get(userleaf.taxon.label)
                summary_list[12] = 'ANI/Placement'
                summary_list[15] = self.aa_percent_msa(
                    msa_dict.get(summary_list[0]))
                summary_list[16] = trans_table_dict.get(summary_list[0])
                fastani_matching_reference = None
                if len(sorted_prefilter_dict) > 0:
                    fastani_matching_reference = sorted_prefilter_dict[0][0]
                    current_ani = all_fastani_dict.get(userleaf.taxon.label).get(
                        fastani_matching_reference).get('ani')
                    current_af = all_fastani_dict.get(userleaf.taxon.label).get(
                        fastani_matching_reference).get('af')

                    taxa_str = ";".join(self.gtdb_taxonomy.get(
                        add_ncbi_prefix(fastani_matching_reference))[:-1])

                    summary_list[1] = self.standardise_taxonomy(taxa_str)
                    summary_list[2] = fastani_matching_reference
                    summary_list[3] = str(
                        self.species_radius.get(fastani_matching_reference))
                    summary_list[4] = ";".join(self.gtdb_taxonomy.get(
                        add_ncbi_prefix(fastani_matching_reference)))
                    current_ani = all_fastani_dict.get(userleaf.taxon.label).get(
                        fastani_matching_reference).get('ani')
                    summary_list[5] = round(current_ani, 2)
                    current_af = all_fastani_dict.get(userleaf.taxon.label).get(
                        fastani_matching_reference).get('af')
                    summary_list[6] = current_af

                    taxa_str = ";".join(self.gtdb_taxonomy.get(
                        add_ncbi_prefix(fastani_matching_reference)))
                    summary_list[1] = self.standardise_taxonomy(
                        taxa_str)

                    summary_list[13] = 'topological placement and ANI have incongruent species assignments'
                    if len(sorted_dict) > 0:
                        other_ref = '; '.join(self._formatnote(
                            sorted_dict, [fastani_matching_reference]))
                        if len(other_ref) == 0:
                            summary_list[14] = None
                        else:
                            summary_list[14] = other_ref

                    summaryfout.write("{}\n".format(
                        '\t'.join(['N/A' if x is None else str(x) for x in summary_list])))

                    classified_user_genomes.append(userleaf.taxon.label)
                else:
                    if len(sorted_dict) > 0:
                        other_ref = '; '.join(self._formatnote(
                            sorted_dict, []))
                        if len(other_ref) == 0:
                            summary_list[14] = None
                        else:
                            summary_list[14] = other_ref
                    unclassified_user_genomes[userleaf.taxon.label] = summary_list
        return classified_user_genomes, unclassified_user_genomes

    def _fastaniWorker(self, sublist_genomes, genomes, out_q):
        """Multi thread worker to calculate FastANI"""
        try:
            for userleaf, potential_nodes in sublist_genomes.iteritems():
                dict_parser_distance = self._calculate_fastani_distance(
                    userleaf, potential_nodes, genomes)
                for k, v in dict_parser_distance.iteritems():
                    if k in out_q:
                        raise Exception("{} not in output.".format(k))
                    out_q[k] = v
            return True
        except Exception as error:
            print(error)
            raise

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
        return ';'.join(subtax)

    def _parse_subnodes(self, list_subnode, closest_rank):
        subtax = []
        multirefrank = []
        initial_loop = True
        for item in list_subnode:
            # We get the taxonomy of all reference genomes
            if item.startswith('RS_') or item.startswith('GB_') or item.startswith('UBA'):
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
        string
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
        for n in tree.preorder_node_iter(lambda n: n != tree.seed_node):
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

    def _calculate_fastani_distance(self, user_leaf, list_leaf, genomes):
        """ Calculate the FastANI distance between all user genomes and the reference to classfy them at the species level

        Parameters
        ----------
        user_leaf : User genome node
        list_leaf : Dictionary of nodes including one or many user genomes and one reference genome.
        genomes : Dictionary of user genomes d[genome_id] -> FASTA file

        Returns
        -------
        dictionary
            dict_results[user_g]={ref_genome1:{"af":af,"ani":ani},ref_genome2:{"af":af,"ani":ani}}
        """
        try:
            self.tmp_output_dir = tempfile.mkdtemp()
            make_sure_path_exists(self.tmp_output_dir)
            dict_parser_distance = {}

            # we first calculate the user genome vs the reference
            # we write the two input files for fastani, the query file and
            # reference file
            path_user_list = os.path.join(
                self.tmp_output_dir, 'query_list.txt')
            with open(path_user_list, 'w') as f:
                f.write('{0}\n'.format(genomes.get(user_leaf.taxon.label)))

            leafnodes = list_leaf.get("potential_g")
            for node in leafnodes:
                leafnode = node[0]
                shortleaf = leafnode.taxon.label
                path_ref_list = os.path.join(self.tmp_output_dir, 'ref_{}.txt'.format(shortleaf))
                if leafnode.taxon.label.startswith('GB_') or leafnode.taxon.label.startswith('RS_'):
                    shortleaf = leafnode.taxon.label[3:]
                with open(path_ref_list, 'w') as f:
                    f.write('{}\n'.format(os.path.join(
                        Config.FASTANI_GENOMES, shortleaf + Config.FASTANI_GENOMES_EXT)))
                # run fastANI
                if not os.path.isfile(path_user_list) or not os.path.isfile(path_ref_list):
                    raise

                path_results = os.path.join(self.tmp_output_dir, 'results_{}_UvsRef.tab'.format(shortleaf))
                path_error = os.path.join(self.tmp_output_dir, 'error_{}.log'.format(shortleaf))

                cmd = 'fastANI --ql {0} --rl {1} -o {2} > /dev/null 2>{3}'.format(path_user_list,
                                                                                  path_ref_list,
                                                                                  path_results,
                                                                                  path_error)

                os.system(cmd)

                if not os.path.isfile(path_results):
                    errstr = 'FastANI has stopped:\n'
                    if os.path.isfile(path_error):
                        with open(path_error) as debug:
                            for line in debug:
                                finalline = line
                            errstr += finalline
                    raise ValueError(errstr)

                dict_parser_distance = self._parse_fastani_results(path_results, dict_parser_distance)

                # We then calculate the reference vs user genome
                path_results_reverse = os.path.join(self.tmp_output_dir, 'results_{}_RefvsU.tab'.format(shortleaf))
                cmd_reverse = 'fastANI --ql {0} --rl {1} -o {2} > /dev/null 2>{3}'.format(path_ref_list,
                                                                                          path_user_list,
                                                                                          path_results_reverse,
                                                                                          path_error)
                os.system(cmd_reverse)
                if not os.path.isfile(path_results_reverse):
                    errstr = 'FastANI has stopped:\n'
                    if os.path.isfile(path_error):
                        with open(path_error) as debug:
                            for line in debug:
                                finalline = line
                            errstr += finalline
                    raise ValueError(errstr)
                dict_parser_distance = self._parse_fastani_results_reverse(path_results_reverse, dict_parser_distance)

            shutil.rmtree(self.tmp_output_dir)
            return dict_parser_distance

        except ValueError as error:
            if os.path.exists(self.tmp_output_dir):
                shutil.rmtree(self.tmp_output_dir)
            raise error
        except Exception as error:
            if os.path.exists(self.tmp_output_dir):
                shutil.rmtree(self.tmp_output_dir)
            raise error

    def _parse_fastani_results(self, fastout_file, dict_results):
        """ Parse the fastani output file


        Parameters
        ----------
        fastout_file : fastani output file.


        Returns
        -------
        dictionary
            dict_results[user_g]={ref_genome1:{"af":af,"ani":ani},ref_genome2:{"af":af,"ani":ani}}
        """
        with open(fastout_file) as fastfile:
            for line in fastfile:
                info = line.strip().split()
                ref_genome = os.path.basename(info[1]).replace(
                    Config.FASTANI_GENOMES_EXT, "")
                user_g = remove_extension(os.path.basename(info[0]))
                ani = float(info[2])
                af = round(float(info[3]) / float(info[4]), 2)
                if user_g in dict_results:
                    dict_results[user_g][ref_genome] = {"ani": ani, 'af': af}
                else:
                    dict_results[user_g] = {ref_genome: {"ani": ani, "af": af}}

        return dict_results

    def _parse_fastani_results_reverse(self, fastout_file, dict_parser_distance):
        # TODO: Merge _parse_fastani_results and _parse_fastani_results_reverse
        """ Parse the fastani output file for the reverse comparison and pick the best ANI and AF


        Parameters
        ----------
        fastout_file : fastani output file.
        dict_parser_distance: dictionaryof user genomes vs list of refrence genomes with ANI and AF


        Returns
        -------
        dictionary
            dict_parser_distance[user_g]={ref_genome1:{"af":af,"ani":ani},ref_genome2:{"af":af,"ani":ani}}
        """
        with open(fastout_file) as fastfile:
            for line in fastfile:
                info = line.strip().split()
                ref_genome = os.path.basename(info[0]).replace(
                    Config.FASTANI_GENOMES_EXT, "")
                user_g = remove_extension(os.path.basename(info[1]))
                ani = float(info[2])
                af = round(float(info[3]) / float(info[4]), 2)
                if user_g in dict_parser_distance:
                    if ref_genome in dict_parser_distance.get(user_g):
                        if dict_parser_distance.get(user_g).get(ref_genome).get('ani') < ani:
                            dict_parser_distance[user_g][ref_genome]["ani"] = ani
                        if dict_parser_distance.get(user_g).get(ref_genome).get('af') < af:
                            dict_parser_distance[user_g][ref_genome]["af"] = af
                    else:
                        dict_parser_distance[user_g][ref_genome] = {"ani": ani, 'af': af}
                else:
                    dict_parser_distance[user_g] = {ref_genome: {"ani": ani, "af": af}}

        return dict_parser_distance

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
        for taxon_id, taxa in taxonomy.iteritems():
            for taxon in taxa:
                taxa_for_dist_inference.add(taxon)

        # sanity check species names as these are a common problem
        species = set()
        for taxon_id, taxa in taxonomy.iteritems():
            if len(taxa) > Taxonomy.rank_index['s__']:
                species_name = taxa[Taxonomy.rank_index['s__']]
                valid, error_msg = True, None
                if species_name != 's__':
                    valid, error_msg = Taxonomy().validate_species_name(
                        species_name, require_full=True, require_prefix=True)
                if not valid:
                    print('[Warning] Species name {} for {} is invalid: {}'.format(species_name, taxon_id, error_msg))
                    continue

                species.add(species_name)

        # restrict taxa to those with a sufficient number of named children
        # Note: a taxonomic group with no children will not end up in the
        # taxon_children data structure so care must be taken when applying
        # this filtering criteria.
        if min_children > 0:
            valid_taxa = set()
            for taxon, children_taxa in taxon_children.iteritems():
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
                    print('[Error] Tree does not contain support values. As such, --min_support should be set to 0.')
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
            sys.stdout.write('{}\r'.format(status_msg))
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
        outgroup : iterable
          Labels of taxa in outgroup.

        Returns
        -------
        Dendropy Tree
            Deep-copy of original tree rerooted on outgroup.
        """

        new_tree = input_tree.clone()

        outgroup = set()
        for genome_id, taxa in taxonomy.iteritems():
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
        # taxa. Finding the MRCA requires a rooted tree and we have no gaurantee
        # that the tree isn't currently rooted within the outgroup clade. There is
        # also no way to identify a node that is gauranteed to be outside the outgroup
        # clade. As such, the tree is randomly rooted on a leaf node not in the outgroup.
        # This random rerooting is performed until the MRCA does not spans all taxa in
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
