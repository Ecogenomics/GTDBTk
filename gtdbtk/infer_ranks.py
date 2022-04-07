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

import dendropy

from gtdbtk.biolib_lite.newick import parse_label, create_label
from gtdbtk.biolib_lite.taxonomy import Taxonomy
from gtdbtk.config.config import (TAXONOMY_FILE,
                                  RED_DIST_BAC_DICT,
                                  RED_DIST_ARC_DICT,
                                  MRCA_RED_BAC120,
                                  MRCA_RED_AR53,
                                  RED_INTERVAL)
from gtdbtk.exceptions import GTDBTkExit
from gtdbtk.relative_distance import RelativeDistance


class InferRanks(object):
    """Establish taxonomic ranks of internal nodes using RED."""

    def __init__(self):
        """Initialize."""

        self.logger = logging.getLogger('timestamp')

    def _get_ingroup_domain(self, ingroup_taxon) -> str:
        """Get domain on ingroup taxon."""

        # read GTDB taxonomy in order to establish domain on ingroup taxon
        gtdb_taxonomy = Taxonomy().read(TAXONOMY_FILE)
        ingroup_domain = None
        for taxa in gtdb_taxonomy.values():
            if ingroup_taxon in taxa:
                ingroup_domain = taxa[Taxonomy.DOMAIN_IDX]

        if ingroup_domain is None:
            raise GTDBTkExit(f'Ingroup taxon {ingroup_taxon} was not found in '
                             f'the GTDB taxonomy.')

        return ingroup_domain

    def _get_median_reds(self, ingroup_domain: str):
        """Get median RED values for domain of ingroup taxon."""

        # get median RED values for domain
        if ingroup_domain == 'd__Bacteria':
            median_reds = RED_DIST_BAC_DICT
        elif ingroup_domain == 'd__Archaea':
            median_reds = RED_DIST_ARC_DICT
        else:
            raise GTDBTkExit(f'Unrecognized GTDB domain: {ingroup_domain}.')

        # report median values
        domain = ingroup_domain.replace('d__', '')
        self.logger.info('Median RED values for {}:'.format(domain))
        for idx, rank_prefix in enumerate(Taxonomy.rank_prefixes):
            if idx != Taxonomy.DOMAIN_IDX and idx != Taxonomy.SPECIES_IDX:
                self.logger.info('  {}\t{:.3f}'.format(
                    Taxonomy.rank_labels[idx].capitalize(),
                    median_reds[rank_prefix]))

        return median_reds

    def _find_ingroup_taxon(self, ingroup_taxon, tree):
        """Find node of ingroup taxon in tree."""

        ingroup_node = None
        for node in tree.postorder_node_iter():
            support, taxon, auxiliary_info = parse_label(node.label)

            if taxon:
                taxa = [t.strip() for t in taxon.split(';')]
                if ingroup_taxon in taxa:
                    if ingroup_node is not None:
                        raise GTDBTkExit(f'Ingroup taxon {ingroup_taxon} '
                                         f'identified multiple times.')
                    ingroup_node = node

        if ingroup_node is None:
            raise GTDBTkExit(f'Ingroup taxon {ingroup_taxon} not found in tree.')

        return ingroup_node

    def _find_ingroup_red(self, ingroup_node, ingroup_domain, tree):
        """Find RED of the ingroup taxon."""

        red_file = MRCA_RED_BAC120
        if ingroup_domain == 'd__Archaea':
            red_file = MRCA_RED_AR53

        # create map from leave labels to tree nodes
        leaf_node_map = {}
        for leaf in tree.leaf_node_iter():
            leaf_node_map[leaf.taxon.label] = leaf

        # find RED value of ingroup node
        reference_nodes = set()
        with open(red_file) as rf:
            for line in rf:
                label_ids, red = line.strip().split('\t')
                labels = label_ids.split('|')
                if len(labels) == 2:
                    taxa = [leaf_node_map[label].taxon for label in labels]
                    node = tree.mrca(taxa=taxa)
                    if node == ingroup_node:
                        return float(red)

        raise GTDBTkExit(f'Could not determine RED of ingroup taxon {ingroup_node}.')

    def _determine_red_ranks(self, node_red, median_reds):
        """Determine suitable taxonomic ranks for node using RED."""

        red_ranks = {}
        for rank_prefix, median_red in median_reds.items():
            rank_idx = Taxonomy.rank_index[rank_prefix]
            rank_label = Taxonomy.rank_labels[rank_idx]

            abs_red_diff = abs(node_red - median_red)
            if abs_red_diff <= RED_INTERVAL:
                red_ranks[rank_label] = abs_red_diff

        red_ranks_label = []
        for rank_label, abs_red_diff in sorted(red_ranks.items(), key=lambda kv: kv[1]):
            red_ranks_label.append(rank_label)

        return '&'.join(red_ranks_label)

    def run(self, input_tree, ingroup_taxon, output_tree):
        """Establish taxonomic ranks of internal nodes using RED..

        Parameters
        ----------
        input_tree : str
          Rooted tree with labelled outgroup.
        ingroup_taxon : str
          Ingroup from which to infer ranks based on RED.
        output_tree: str
          Output directory.
        """

        # get domain on ingroup taxon
        ingroup_domain = self._get_ingroup_domain(ingroup_taxon)

        # get median RED values for domain of ingroup taxon
        median_reds = self._get_median_reds(ingroup_domain)

        # read tree
        self.logger.info('Reading tree.')
        tree = dendropy.Tree.get_from_path(input_tree,
                                           schema='newick',
                                           rooting='force-rooted',
                                           preserve_underscores=True)

        # find ingroup taxon
        ingroup_node = self._find_ingroup_taxon(ingroup_taxon, tree)

        # get RED of ingroup taxon
        ingroup_red = self._find_ingroup_red(ingroup_node, ingroup_domain, tree)
        self.logger.info('RED of ingroup taxon {} = {:.3f}'.format(
            ingroup_taxon, ingroup_red))

        # get RED value of ingroup taxon
        self.logger.info('Decorating tree with RED and rank information.')
        red = RelativeDistance()
        red.decorate_rel_dist(ingroup_node, ingroup_red)

        for node in ingroup_node.preorder_iter():
            if node.is_leaf():
                continue

            support, taxon, auxiliary_info = parse_label(node.label)

            if auxiliary_info:
                auxiliary_info += '|RED={:.3f}'.format(node.rel_dist)
            else:
                auxiliary_info = 'RED={:.3f}'.format(node.rel_dist)

            red_ranks = self._determine_red_ranks(node.rel_dist, median_reds)
            auxiliary_info += '|{}'.format(red_ranks)

            new_label = create_label(support, taxon, auxiliary_info)
            node.label = new_label

        # write RED decorated tree to file
        tree.write_to_path(output_tree,
                           schema='newick',
                           suppress_rooting=True,
                           unquoted_underscores=True)
