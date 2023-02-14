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

import gtdbtk.config.config as Config
from gtdbtk.biolib_lite.common import make_sure_path_exists
from gtdbtk.biolib_lite.newick import parse_label
from gtdbtk.biolib_lite.seq_io import read_fasta
from gtdbtk.config.output import *
from gtdbtk.exceptions import GenomeMarkerSetUnknown, GTDBTkExit
from gtdbtk.files.classify_summary import ClassifySummaryFileRow
from gtdbtk.files.pplacer_classification import PplacerHighClassifyRow, PplacerHighClassifyFile
from gtdbtk.files.tree_mapping import GenomeMappingFileRow
from gtdbtk.tools import TreeTraversal, standardise_taxonomy, aa_percent_msa

class Split(object):
    """Determine taxonomic classification of genomes by ML placement using the Split Methods."""

    def __init__(self, order_rank, gtdb_taxonomy, reference_ids):
        """Initialize."""
        self.order_rank = order_rank
        self.logger = logging.getLogger('timestamp')
        self.gtdb_taxonomy = gtdb_taxonomy
        self.reference_ids = reference_ids

        # rank_of_interest determine the rank in the tree_mapping file for
        # lower classification
        self.rank_of_interest = "c__"

    def get_high_pplacer_taxonomy(self, out_dir, marker_set_id, prefix, user_msa_file, tree):
        """Parse the pplacer tree and write the partial taxonomy for each user genome based on their placements

        Parameters
        ----------
        out_dir : output directory
        prefix : desired prefix for output files
        marker_set_id : bacterial or archaeal id (bac120 or ar53)
        user_msa_file : msa file listing all user genomes for a certain domain
        tree : pplacer tree including the user genomes

        Returns
        -------
        dictionary[genome_label]=pplacer_taxonomy

        """
        results = {}
        out_root = os.path.join(out_dir, 'classify', 'intermediate_results')
        make_sure_path_exists(out_root)

        if marker_set_id == 'bac120':
            out_pplacer = PplacerHighClassifyFile(out_dir,prefix)
        else:
            self.logger.error('There was an error determining the marker set.')
            raise GenomeMarkerSetUnknown

        red_bac_dict = Config.RED_DIST_BAC_DICT

        # We get the pplacer taxonomy for comparison
        user_genome_ids = set(read_fasta(user_msa_file).keys())
        for leaf in tree.leaf_node_iter():

            is_on_terminal_branch = False
            terminal_branch_test = False
            term_branch_taxonomy = ''
            if leaf.taxon.label in user_genome_ids:
                pplacer_row = PplacerHighClassifyRow()
                taxa = []
                cur_node = leaf
                current_rel_dist = 1.0
                # every user genomes has a RED value of one assigned to it
                while cur_node.parent_node:
                    # we go up the tree from the user genome
                    if hasattr(cur_node, 'rel_dist') and current_rel_dist == 1.0 and cur_node.rel_dist < 1.0:
                        # if the parent node of the current genome has a red distance,
                        # it means it is part of the reference tree
                        # we store the first RED value encountered in the
                        # tree
                        current_rel_dist = cur_node.rel_dist
                    if cur_node.is_internal():
                        # We check if the genome is place on a terminal
                        # branch

                        if not terminal_branch_test:
                            child_genomes = [nd.taxon.label for nd in cur_node.leaf_nodes(
                            ) if nd.taxon.label not in user_genome_ids]
                            if len(child_genomes) == 1:
                                is_on_terminal_branch = True
                                term_branch_taxonomy = self.gtdb_taxonomy.get(
                                    child_genomes[0])
                                terminal_branch_test = True
                            if len(child_genomes) > 1:
                                terminal_branch_test = True
                    # While going up the tree we store of taxonomy
                    # information
                    _support, taxon, _aux_info = parse_label(
                        cur_node.label)
                    if taxon:
                        for t in taxon.split(';')[::-1]:
                            taxa.append(t.strip())
                    cur_node = cur_node.parent_node

                taxa_str = ';'.join(taxa[::-1])

                pplacer_tax = str(taxa_str)

                taxa_str_terminal,taxa_str_red = '',''

                if is_on_terminal_branch:
                    # some rank may be missing from going up the tree.
                    # if the genome is on a terminal branch,
                    # we can select the taxonomy from the reference leaf to get the low level of the taxonomy
                    # we select down to genus
                    if len(taxa) > 1:
                        tax_of_leaf = term_branch_taxonomy[term_branch_taxonomy.index(
                            taxa_str.split(';')[-1]) + 1:-1]
                    else:
                        tax_of_leaf = term_branch_taxonomy[1:-1]
                        taxa_str = 'd__Bacteria'

                    taxa_str_terminal = self._classify_on_terminal_branch(
                        tax_of_leaf, current_rel_dist, taxa_str.split(';')[-1][0:3], term_branch_taxonomy,
                        red_bac_dict)

                cur_node = leaf
                parent_taxon_node = cur_node.parent_node
                _support, parent_taxon, _aux_info = parse_label(
                    parent_taxon_node.label)

                while parent_taxon_node is not None and not parent_taxon:
                    parent_taxon_node = parent_taxon_node.parent_node
                    _support, parent_taxon, _aux_info = parse_label(
                        parent_taxon_node.label)

                # is the node represent multiple ranks, we select the lowest one
                # i.e. if node is p__A;c__B;o__C we pick o__
                parent_rank = parent_taxon.split(";")[-1]

                if parent_rank[0:3] != 'g__':
                    node_in_ref_tree = cur_node
                    while len([childnd.taxon.label.replace("'", '') for childnd in node_in_ref_tree.leaf_iter(
                    ) if childnd.taxon.label in self.reference_ids]) == 0:
                        node_in_ref_tree = node_in_ref_tree.parent_node
                    # we select a node of the reference tree

                    # we select the child rank (if parent_rank = 'c__'
                    # child rank will be 'o__)'
                    child_rk = self.order_rank[self.order_rank.index(
                        parent_rank[0:3]) + 1]

                    # get all reference genomes under the current node
                    list_subnode = [childnd.taxon.label.replace("'", '') for childnd in
                                    node_in_ref_tree.leaf_iter()
                                    if childnd.taxon.label in self.reference_ids]

                    # get all names for the child rank
                    list_ranks = [self.gtdb_taxonomy.get(name)[self.order_rank.index(child_rk)]
                                  for name in list_subnode]

                    # if there is just one rank name
                    if len(set(list_ranks)) == 1:
                        child_taxons = []
                        child_rel_dist = None
                        for subranknd in node_in_ref_tree.preorder_iter():
                            _support, subranknd_taxon, _aux_info = parse_label(
                                subranknd.label)
                            if subranknd.is_internal() and subranknd_taxon is not None and subranknd_taxon.startswith(
                                    child_rk):
                                child_taxons = subranknd_taxon.split(
                                    ";")
                                child_taxon_node = subranknd
                                child_rel_dist = child_taxon_node.rel_dist
                                break

                        taxa_str_red, taxa_str_terminal = self._classify_on_internal_branch(leaf.taxon.label,
                                                                                            child_taxons,
                                                                                            current_rel_dist,
                                                                                            child_rel_dist,
                                                                                            node_in_ref_tree,
                                                                                            parent_rank, child_rk,
                                                                                            taxa_str,
                                                                                            taxa_str_terminal,
                                                                                            is_on_terminal_branch,
                                                                                            red_bac_dict)
                    else:
                        taxa_str_red = taxa_str


                results[leaf.taxon.label] = {"tk_tax_red": standardise_taxonomy(taxa_str_red, 'bac120'),
                                             "tk_tax_terminal": standardise_taxonomy(taxa_str_terminal,
                                                                                     'bac120'),
                                             "pplacer_tax": standardise_taxonomy(pplacer_tax, 'bac120'),
                                             'rel_dist': current_rel_dist}

                pplacer_row.gid = leaf.taxon.label
                pplacer_row.gtdb_taxonomy_red = standardise_taxonomy(taxa_str_red, 'bac120')
                pplacer_row.gtdb_taxonomy_terminal = standardise_taxonomy(taxa_str_terminal, 'bac120')
                pplacer_row.pplacer_taxonomy = standardise_taxonomy(pplacer_tax, 'bac120')
                pplacer_row.is_terminal = is_on_terminal_branch
                pplacer_row.red = current_rel_dist

                out_pplacer.add_row(pplacer_row)

        out_pplacer.write()
        return results

    def _classify_on_internal_branch(self, leaf, child_taxons, current_rel_list, child_rel_dist, node_in_ref_tree,
                                     parent_rank, child_rk, taxa_str, taxa_str_terminal, is_on_terminal_branch,
                                     red_bac_dict):
        """
         Classification on an internal node is very similar to the 'normal' classification
         """

        # Persist descendant information for efficient traversal.
        tt = TreeTraversal()

        closest_rank = None

        if len(child_taxons) == 0:
            list_leaves = [childnd.taxon.label.replace("'", '')
                           for childnd in tt.get_leaf_nodes(node_in_ref_tree)
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
                    raise GTDBTkExit('There should be only one leaf.')
                else:
                    closest_rank = parent_rank
                    detection = "taxonomic classification fully defined by topology"
            list_leaf_ranks = self.gtdb_taxonomy.get(list_leaves[0])[
                              self.order_rank.index(child_rk):-1]  # We remove the species name

            for leaf_taxon in reversed(list_leaf_ranks):
                leaf_taxon_rank = leaf_taxon[:3]
                if leaf_taxon == list_leaf_ranks[0]:
                    if abs(current_rel_list - red_bac_dict.get(leaf_taxon_rank)) < abs(
                            current_rel_list - red_bac_dict.get(parent_rank[:3])):
                        closest_rank = leaf_taxon
                        break
                else:
                    pchildrank = list_leaf_ranks[list_leaf_ranks.index(leaf_taxon) - 1]
                    if abs(current_rel_list - red_bac_dict.get(leaf_taxon_rank)) < abs(
                            current_rel_list - red_bac_dict.get(pchildrank[:3])):
                        closest_rank = leaf_taxon
                        break
            if closest_rank is None:
                closest_rank = parent_rank
        # if there is multiple ranks on the child node (i.e genome between p__Nitrospirae and c__Nitrospiria;o__Nitrospirales;f__Nitropiraceae)
        # we loop through the list of rank from f_ to c_ rank
        for child_taxon in reversed(child_taxons):
            child_taxon_rank = child_taxon[:3]
            if child_taxon == child_taxons[0]:
                if (abs(current_rel_list - red_bac_dict.get(child_taxon_rank)) < abs(
                        child_rel_dist - red_bac_dict.get(child_taxon_rank)) and
                        abs(current_rel_list - red_bac_dict.get(child_taxon_rank)) < abs(
                            current_rel_list - red_bac_dict.get(parent_rank[:3]))):
                    closest_rank = child_taxon
                elif closest_rank is None:
                    closest_rank = parent_rank
            else:
                pchildrank = child_taxons[child_taxons.index(
                    child_taxon) - 1]
                if (abs(current_rel_list - red_bac_dict.get(child_taxon_rank)) < abs(
                        current_rel_list - red_bac_dict.get(child_taxon_rank)) and
                        abs(current_rel_list - red_bac_dict.get(child_taxon_rank)) < abs(
                            child_rel_dist - red_bac_dict.get(child_taxon_rank))):
                    closest_rank = child_taxon
                    break
        if closest_rank is not None:
            # when we have the closest rank found, we can find it in
            # gtdb_Taxonomy and get the higher level from it.
            for k, v in self.gtdb_taxonomy.items():
                if closest_rank in v:
                    taxa_str = ';'.join(v[1:v.index(closest_rank) + 1])
                    # All classification should be at least to the class level if a genome
                    # is placed on an internal branch with only one class under
                    if any(x.startswith('c__') for x in child_taxons) \
                            and self.order_rank.index(closest_rank[0:3]) < self.order_rank.index('c__') \
                            and ('c__' in taxa_str_terminal.split(';') or not is_on_terminal_branch):
                        taxa_str_terminal = ';'.join(v[1:self.order_rank.index('c__') + 1])
                    break

        return taxa_str, taxa_str_terminal

    def _classify_on_terminal_branch(self, list_leaf_ranks, current_rel_list, parent_rank, term_branch_taxonomy,
                                     red_bac_dict):
        """
        When a genome is on a terminal branch, we can guess the low level of its taxonomy,
        based on the RED value
        :param list_leaf_ranks: Taxonomy of the reference leaf
        :param current_rel_list: RED value for the genome of interest
        :param parent_rank: Parent  rank of the genome of interest
        :param term_branch_taxonomy: Full taxonomy of the branch , down to genus
        :param red_bac_dict: RED dictionary for Bacteria
        :return: taxonomy for the genome of interest based on RED
        """
        closest_rank = None
        for leaf_taxon in reversed(list_leaf_ranks):
            if leaf_taxon == list_leaf_ranks[0]:
                if abs(current_rel_list - red_bac_dict.get(leaf_taxon[:3])) < abs(
                        current_rel_list - red_bac_dict.get(parent_rank)):
                    closest_rank = leaf_taxon[:3]
                    break
            else:
                pchildrank = list_leaf_ranks[list_leaf_ranks.index(
                    leaf_taxon) - 1]
                if abs(current_rel_list - red_bac_dict.get(leaf_taxon[:3])) < abs(
                        current_rel_list - red_bac_dict.get(pchildrank[:3])):
                    closest_rank = leaf_taxon[:3]
                    break
        if closest_rank is None:
            closest_rank = parent_rank
        # All classification should be at least to the class level if a genome
        # is placed on a terminal branch
        if self.order_rank.index(closest_rank) < self.order_rank.index('c__'):
            return ';'.join(term_branch_taxonomy[1:self.order_rank.index('c__') + 1])

        return ';'.join(term_branch_taxonomy[1:self.order_rank.index(closest_rank) + 1])

    def map_high_taxonomy(self,high_classification, mapping_dict, summary_file,
                          tree_mapping_file,msa_dict,trans_table_dict,percent_multihit_dict,bac_ar_diff,warning_counter):
        mapped_rank = {}
        counter = 0
        high_taxonomy_used = {}
        for k, v in high_classification.items():
            # if the classification has an order
            rk_to_check = v.get('tk_tax_red').split(
                ';')[self.order_rank.index(self.rank_of_interest)]
            if len(rk_to_check) > 3:
                mapped_rank.setdefault(
                    mapping_dict.get(rk_to_check), []).append(k)
                counter += 1
                high_taxonomy_used[k] = ["RED",v.get('tk_tax_terminal'),v.get('tk_tax_red')]
            else:
                rk_to_check = v.get('tk_tax_terminal').split(
                    ';')[self.order_rank.index(self.rank_of_interest)]
                if len(rk_to_check) > 3:
                    mapped_rank.setdefault(
                        mapping_dict.get(rk_to_check), []).append(k)
                    counter += 1
                    high_taxonomy_used[k] = ["TERMINAL",v.get('tk_tax_terminal'),v.get('tk_tax_red')]
                else:
                    summary_row = ClassifySummaryFileRow()
                    summary_row.gid = k
                    summary_row.classification = v.get('tk_tax_red')
                    summary_row.pplacer_tax = v.get('pplacer_tax')
                    summary_row.red_value = v.get('rel_dist')
                    summary_row.note = 'classification based on placement in backbone tree'
                    summary_row.msa_percent = aa_percent_msa(msa_dict.get(summary_row.gid))
                    summary_row.tln_table = trans_table_dict.get(summary_row.gid)

                    warnings = []
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


                    mapping_row = GenomeMappingFileRow()
                    mapping_row.gid = k
                    mapping_row.ani_classification = False
                    mapping_row.mapped_tree = 'backbone'
                    mapping_row.rule = 'Rule 1'

                    tree_mapping_file.add_row(mapping_row)
                    summary_file.add_row(summary_row)

        return mapped_rank,warning_counter, counter , high_taxonomy_used
