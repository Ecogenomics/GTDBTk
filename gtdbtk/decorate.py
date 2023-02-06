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
import sys
from collections import defaultdict, namedtuple

import dendropy

from gtdbtk.biolib_lite.newick import parse_label, create_label
from gtdbtk.biolib_lite.taxonomy import Taxonomy


class Decorate(object):
    """Decorate internal nodes with taxa labels."""

    def __init__(self):
        """Initialize."""

        self.logger = logging.getLogger('timestamp')

        self.StatsTable = namedtuple('StatsTable',
                                     'node fmeasure precision recall taxa_in_lineage total_taxa num_leaves_with_taxa rogue_out rogue_in')

    def _fmeasure(self, tree, taxonomy):
        """Find node with highest F-measure for each taxon.
        
        Finds best placement for each taxon label
        by calculating the F-measure for every taxon
        at every node.
        
        Parameters
        ----------
        tree : Tree
          Dendropy Tree.
        taxonomy : d[extent_taxon_id] -> taxa list
          Taxon labels for extant taxa.
          
        Returns
        -------
        d[taxon] -> [(Node, F-measure, precision, recall_, ...]
            Node(s) with highest F-measure for each taxon.
        """

        # get named lineages/taxa at each taxonomic rank
        taxa_at_rank = Taxonomy().named_lineages_at_rank(taxonomy)

        # get extant taxa for each taxon label
        extent_taxa_with_label = {}
        for i, rank in enumerate(Taxonomy.rank_labels):
            extent_taxa_with_label[i] = Taxonomy().extant_taxa_for_rank(rank, taxonomy)

        # get parent taxon for each taxon:
        taxon_parents = Taxonomy().parents(taxonomy)

        # get number of leaves and taxon in each lineage
        self.logger.info('Calculating taxa within each lineage.')
        for node in tree.preorder_node_iter():
            num_leaves = 0
            taxa_count = defaultdict(lambda: defaultdict(int))
            for leaf in node.leaf_iter():
                num_leaves += 1
                for rank_index, taxon in enumerate(taxonomy[leaf.taxon.label]):
                    if taxon != Taxonomy.rank_prefixes[rank_index]:
                        taxa_count[rank_index][taxon] += 1

            node.num_leaves = num_leaves
            node.taxa_count = taxa_count

        taxa_in_tree = defaultdict(int)
        for leaf in tree.leaf_node_iter():
            for taxon in taxonomy[leaf.taxon.label]:
                taxa_in_tree[taxon] += 1

        # find node with best F-measure for each taxon
        fmeasure_for_taxa = {}
        for rank_index in range(0, len(Taxonomy.rank_labels)):
            # if rank_index == 6: #*** skip species
            #    continue 
            self.logger.info('Processing {:,} taxa at {} rank.'.format(
                len(taxa_at_rank[rank_index]),
                Taxonomy.rank_labels[rank_index].capitalize()))

            for taxon in taxa_at_rank[rank_index]:
                if rank_index == 0:
                    # processing taxa at the domain is a special case
                    taxon_parent_node = tree.seed_node
                else:
                    # find first named parent 
                    # e.g., Cyanobacteria for Synechococcales in d__Bacteria;p__Cyanobacteria;c__;o__Synechococcales
                    parent_taxon = 'x__'
                    parent_index = rank_index - 1
                    while len(parent_taxon) == 3 and parent_index != -1:
                        parent_taxon = taxon_parents[taxon][parent_index]
                        parent_index -= 1

                    if parent_taxon in fmeasure_for_taxa:
                        # only need to process the lineage below the parent node,
                        # but must take the MRCA if the placement of the parent
                        # taxon is unresolved
                        parent_nodes = []
                        for stat_table in fmeasure_for_taxa[parent_taxon]:
                            parent_nodes.append(stat_table.node)

                        if len(parent_nodes) == 1:
                            taxon_parent_node = parent_nodes[0]
                        else:
                            taxa = []
                            for p in parent_nodes:
                                taxa += [leaf.taxon for leaf in p.leaf_iter()]
                            taxon_parent_node = tree.mrca(taxa=taxa)

                        if taxon_parent_node.taxa_count[rank_index][taxon] < 0.5 * taxa_in_tree[taxon]:
                            # substantial portion of genomes for taxon fall outside 
                            # the parent lineages so best search the entire tree
                            taxon_parent_node = tree.seed_node
                    else:
                        # the parent for this taxon was not placed so
                        # it can be ignored (e.g., bacterial phylum in archaeal tree)
                        continue

                cur_taxon_fmeasure = -1
                cur_taxa = set(extent_taxa_with_label[rank_index][taxon])
                total_taxa = len(cur_taxa)

                for node in taxon_parent_node.preorder_iter():
                    taxa_in_lineage = node.taxa_count[rank_index][taxon]
                    num_leaves_with_taxa = sum(node.taxa_count[rank_index].values())

                    if taxa_in_lineage != 0 and num_leaves_with_taxa != 0:
                        precision = float(taxa_in_lineage) / num_leaves_with_taxa
                        recall = float(taxa_in_lineage) / total_taxa
                        fmeasure = (2 * precision * recall) / (precision + recall)

                        if fmeasure >= cur_taxon_fmeasure:
                            node_taxa = set([l.taxon.label for l in node.leaf_iter()])
                            rogue_out = cur_taxa - node_taxa
                            rogue_in = []
                            for gid in node_taxa - cur_taxa:
                                if taxonomy[gid][rank_index] != Taxonomy.rank_prefixes[rank_index]:
                                    rogue_in.append(gid)

                            stat_table = self.StatsTable(node=node,
                                                         fmeasure=fmeasure,
                                                         precision=precision,
                                                         recall=recall,
                                                         taxa_in_lineage=taxa_in_lineage,
                                                         total_taxa=total_taxa,
                                                         num_leaves_with_taxa=num_leaves_with_taxa,
                                                         rogue_out=rogue_out,
                                                         rogue_in=rogue_in)

                            if fmeasure > cur_taxon_fmeasure:
                                cur_taxon_fmeasure = fmeasure
                                fmeasure_for_taxa[taxon] = [stat_table]
                            elif fmeasure == cur_taxon_fmeasure:
                                fmeasure_for_taxa[taxon].append(stat_table)

        return fmeasure_for_taxa

    def _strip_taxon_labels(self, tree):
        """Remove any previous taxon labels.
        
        Parameters
        ----------
        tree : Tree
          Dendropy Tree.
        """

        for node in tree.internal_nodes():
            support, _taxon, _aux_info = parse_label(node.label)
            if support is not None:
                node.label = create_label(support, None, None)
            else:
                node.label = None

    def _assign_taxon_labels(self, fmeasure_for_taxa):
        """Assign taxon labels to nodes.
        
        Parameters
        ----------
        fmeasure_for_taxa : d[taxon] -> [(Node, F-measure, precision, recall), ...]
          Node with highest F-measure for each taxon.
          
        Returns
        -------
        set
            Taxon labels placed in tree.
        """

        placed_taxon = set()
        for taxon in Taxonomy().sort_taxa(fmeasure_for_taxa.keys()):
            if len(fmeasure_for_taxa[taxon]) == 1:
                placed_taxon.add(taxon)

                stat_table = fmeasure_for_taxa[taxon][0]
                node = stat_table.node
                fmeasure = stat_table.fmeasure
                precision = stat_table.precision
                recall = stat_table.recall

                support, taxon_label, aux_info = parse_label(node.label)
                if taxon_label:
                    taxon_label += '; ' + taxon
                else:
                    taxon_label = taxon
                node.label = create_label(support, taxon_label, aux_info)

        return placed_taxon

    def _write_statistics_table(self, fmeasure_for_taxa, taxonomy, out_table):
        """Write table containing statistics for each taxon.
        
        Parameters
        ----------
        fmeasure_for_taxa : d[taxon] -> [(Node, F-measure, precision, recall)]
          Node with highest F-measure for each taxon.
        taxonomy : d[unique_id] -> [d__<taxon>; ...; s__<taxon>]
          Taxonomic information for taxa in tree of interest.
        out_table : str
          Output table to write statistics for assigned labels.  
        """

        # get extent taxa
        extant_taxa = Taxonomy().extant_taxa(taxonomy)

        fout_table = open(out_table, 'w')
        fout_table.write('Taxon\tNo. Expected in Tree\tF-measure\tPrecision\tRecall')
        fout_table.write('\tNo. Genomes from Taxon\tNo. Genome In Lineage')
        fout_table.write('\tRogue out\tRogue in\n')
        for taxon in Taxonomy().sort_taxa(fmeasure_for_taxa.keys()):
            if len(fmeasure_for_taxa[taxon]) != 1:
                self.logger.error('Multiple positions specified for taxon label.')
                sys.exit()

            num_genomes = len(extant_taxa[taxon])

            stat_table = fmeasure_for_taxa[taxon][0]
            fout_table.write('%s\t%d\t%.4f\t%.4f\t%.4f\t%d\t%d\t%s\t%s\n' % (
                taxon,
                num_genomes,
                stat_table.fmeasure,
                stat_table.precision,
                stat_table.recall,
                stat_table.taxa_in_lineage,
                stat_table.num_leaves_with_taxa,
                ','.join(stat_table.rogue_out),
                ','.join(stat_table.rogue_in)))

        fout_table.close()

    def _leaf_taxa(self, leaf):
        """Get taxonomic information for leaf node.
        
        Parameters
        ----------
        leaf : Node
          Node in tree.
          
        Returns
        -------
        list
          Taxa for leaf in rank order.
        """

        leaf_taxa = []

        parent = leaf
        while parent:
            _support, taxon, _aux_info = parse_label(parent.label)

            if taxon:
                for t in taxon.split(';')[::-1]:
                    leaf_taxa.append(t.strip())

            parent = parent.parent_node

        ordered_taxa = leaf_taxa[::-1]

        # fill in missing ranks
        last_rank = ordered_taxa[-1][0:3]
        for i in range(Taxonomy.rank_prefixes.index(last_rank) + 1, len(Taxonomy.rank_prefixes)):
            ordered_taxa.append(Taxonomy.rank_prefixes[i])

        return ordered_taxa

    def _write_taxonomy(self, tree, out_taxonomy):
        """Write taxonomy decorated on tree to file.
        
        Parameters
        ----------
        tree : Tree
          Dendropy Tree.
        out_taxonomy : str
          Output file.
        """

        fout = open(out_taxonomy, 'w')
        for leaf in tree.leaf_node_iter():
            taxa = self._leaf_taxa(leaf)
            fout.write('%s\t%s\n' % (leaf.taxon.label, '; '.join(taxa)))

        fout.close()

    def run(self, input_tree, input_taxonomy, output_tree):
        """Decorate internal nodes with taxa labels.

        Parameters
        ----------
        input_tree : str
          Tree to decorate
        input_taxonomy : dict
          Dictionary indicating taxa for genomes.
        output_tree: str
          Name of output tree.
        """

        # read tree
        self.logger.info('Reading tree.')
        tree = dendropy.Tree.get_from_path(input_tree,
                                           schema='newick',
                                           rooting='force-rooted',
                                           preserve_underscores=True)

        # remove any previous taxon labels
        self.logger.info('Removing any previous internal node labels.')
        self._strip_taxon_labels(tree)

        # trim taxonomy to taxa in tree
        taxonomy = {}
        for leaf in tree.leaf_node_iter():
            taxonomy[leaf.taxon.label] = input_taxonomy.get(leaf.taxon.label, Taxonomy.rank_prefixes)

        # find best placement for each taxon based on the F-measure statistic
        self.logger.info('Calculating F-measure statistic for each taxa.')
        fmeasure_for_taxa = self._fmeasure(tree, taxonomy)

        # select most terminal placement for each taxon in order to be conservative
        ambiguous_placements = set()
        for taxon, fmeasures in fmeasure_for_taxa.items():
            if len(fmeasures) != 1:
                ambiguous_placements.add(taxon)
                fmeasure_for_taxa[taxon] = [fmeasures[-1]]

        if len(ambiguous_placements) > 0:
            self.logger.warning('There are {:,} taxa with multiple placements '
                                'of equal quality.'.format(len(ambiguous_placements)))
            self.logger.warning('These were resolved by placing the label at '
                                'the most terminal position.')
            self.logger.warning('Ideally, taxonomic assignment of all genomes '
                                'should be established before tree decoration.')

        # place all labels on tree
        self.logger.info('Placing labels on tree.')
        placed_taxon = self._assign_taxon_labels(fmeasure_for_taxa)

        # write statistics for placed taxon labels
        self.logger.info('Writing out statistics for taxa.')
        out_table = output_tree + '-table'
        self._write_statistics_table(fmeasure_for_taxa, taxonomy, out_table)

        # output taxonomy of extant taxa on tree
        self.logger.info('Writing out inferred taxonomy for each genome.')
        out_taxonomy = output_tree + '-taxonomy'
        self._write_taxonomy(tree, out_taxonomy)

        # output decorated tree
        self.logger.info('Writing out decorated tree.')
        tree.write_to_path(output_tree,
                           schema='newick',
                           suppress_rooting=True,
                           unquoted_underscores=True)

        return {'all':[out_table,out_taxonomy,output_tree]}
