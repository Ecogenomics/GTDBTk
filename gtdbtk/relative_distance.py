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

from collections import defaultdict

import dendropy
from numpy import (median as np_median,
                   array as np_array,
                   arange as np_arange,
                   percentile as np_percentile,
                   ones_like as np_ones_like,
                   histogram as np_histogram)

from gtdbtk.biolib_lite.newick import parse_label
from gtdbtk.biolib_lite.taxonomy import Taxonomy


class RelativeDistance(object):
    """Determine relative rates of evolutionary divergence.

    This code is based on Phylorank: https://github.com/dparks1134/PhyloRank 

    """

    def __init__(self):
        """Initialization."""

        pass

    def _avg_descendant_rate(self, root_node):
        """Calculate average rate of divergence for each nodes in a tree.

        The average rate is the arithmetic mean of the
        branch length to all descendant taxa.

        Parameters
        ----------
        root_node : Dendropy Node
            Root node defining tree or subtree.

        Returns
        -------
        The following attributes are added to each node:
          mean_dist: mean distance to tips
          num_taxa: number of terminal taxa
        """

        # calculate the mean branch length to extant taxa
        for node in root_node.postorder_iter():
            avg_div = 0
            if node.is_leaf():
                node.mean_dist = 0.0
                node.num_taxa = 1
            else:
                node.num_taxa = sum([1 for _ in node.leaf_iter()])
                for c in node.child_node_iter():
                    num_tips = c.num_taxa
                    avg_div += (float(c.num_taxa) / node.num_taxa) * \
                               (c.mean_dist + c.edge_length)

            node.mean_dist = avg_div

    def decorate_rel_dist(self, root_node, root_red=0.0):
        """Calculate relative distance to each internal node.

        Parameters
        ----------
        root_node : Dendropy Node
            Root node defining tree or subtree.

        Returns
        -------
        The following attributes are added to each node:
          mean_dist: mean distance to tips
          num_taxa: number of terminal taxa
          rel_dists: relative distance of node between root and extant organisms
        """

        if isinstance(root_node, dendropy.Tree):
            root_node = root_node.seed_node

        self._avg_descendant_rate(root_node)

        for node in root_node.preorder_iter():
            if node == root_node:
                node.rel_dist = root_red
            elif node.is_leaf():
                node.rel_dist = 1.0
            else:
                a = node.edge_length
                b = node.mean_dist
                x = node.parent_node.rel_dist

                if (a + b) != 0:
                    rel_dist = x + (a / (a + b)) * (1.0 - x)
                else:
                    # internal node has zero length to parent,
                    # so should have the same relative distance
                    # as the parent node
                    rel_dist = x

                node.rel_dist = rel_dist

    def rel_dist_to_named_clades(self, tree):
        """Determine relative distance to specific taxa.

        Parameters
        ----------
        tree : Dendropy Tree
            Phylogenetic tree.

        Returns
        -------
        dict : d[rank_index][taxon] -> relative divergence
        """

        # calculate relative distance for all nodes
        self.decorate_rel_dist(tree)

        # assign internal nodes with ranks from
        rel_dists = defaultdict(dict)
        for node in tree.preorder_node_iter(lambda n: n != tree.seed_node):
            if not node.label or node.is_leaf():
                continue

            # check for support value
            _support, taxon_name, _auxiliary_info = parse_label(node.label)

            if not taxon_name:
                continue

            # get most-specific rank if a node represents multiple ranks
            if ';' in taxon_name:
                taxon_name = taxon_name.split(';')[-1].strip()

            most_specific_rank = taxon_name[0:3]
            rel_dists[Taxonomy.rank_index[most_specific_rank]
                      ][taxon_name] = node.rel_dist

        return rel_dists

    def _distribution_summary_plot(self, phylum_rel_dists, taxa_for_dist_inference, plot_file):
        """Summary plot showing the distribution of taxa at each taxonomic rank under different rootings.

        Parameters
        ----------
        phylum_rel_dists: phylum_rel_dists[phylum][rank_index][taxon] -> relative divergences
            Relative divergence of taxon at each rank for different phylum-level rootings.
        taxa_for_dist_inference : iterable
            Taxa to considered when inferring distributions.
        plot_file : str
            Desired name of output plot.
        """

        self.fig.clear()
        self.fig.set_size_inches(12, 6)
        ax = self.fig.add_subplot(111)

        # determine median relative distance for each taxa
        medians_for_taxa = self.taxa_median_rd(phylum_rel_dists)

        # create percentile and classification boundary lines
        percentiles = {}
        for i, rank in enumerate(sorted(medians_for_taxa.keys())):
            v = [np_median(dists) for taxon, dists in medians_for_taxa[rank].items(
            ) if taxon in taxa_for_dist_inference]
            if not v:
                # not taxa at rank suitable for creating classification
                # boundaries
                continue

            p10, p50, p90 = np_percentile(v, [10, 50, 90])
            ax.plot((p10, p10), (i, i + 0.25),
                    c=(0.3, 0.3, 0.3), lw=2, zorder=2)
            ax.plot((p50, p50), (i, i + 0.5),
                    c=(0.3, 0.3, 0.3), lw=2, zorder=2)
            ax.plot((p90, p90), (i, i + 0.25),
                    c=(0.3, 0.3, 0.3), lw=2, zorder=2)

            for b in [-0.2, -0.1, 0.1, 0.2]:
                boundary = p50 + b
                if 1.0 > boundary > 0.0:
                    if abs(b) == 0.1:
                        c = (1.0, 0.65, 0.0)  # orange
                    else:
                        c = (1.0, 0.0, 0.0)
                    ax.plot((boundary, boundary),
                            (i, i + 0.5), c=c, lw=2, zorder=2)

            percentiles[i] = [p10, p50, p90]

        # create scatter plot and results table
        x = []
        y = []
        c = []
        labels = []
        rank_labels = []
        for i, rank in enumerate(sorted(medians_for_taxa.keys())):
            rank_label = Taxonomy.rank_labels[rank]
            rank_labels.append(rank_label + ' (%d)' %
                               len(medians_for_taxa[rank]))

            mono = []
            poly = []
            no_inference = []
            for clade_label, dists in medians_for_taxa[rank].items():
                md = np_median(dists)
                x.append(md)
                y.append(i)
                labels.append(clade_label)

                if self._is_integer(clade_label.split('^')[-1]):
                    # taxa with a numerical suffix after a caret indicate
                    # polyphyletic groups when decorated with tax2tree
                    c.append((1.0, 0.0, 0.0))
                    poly.append(md)
                elif clade_label not in taxa_for_dist_inference:
                    c.append((0.3, 0.3, 0.3))
                    no_inference.append(md)
                else:
                    c.append((0.0, 0.0, 1.0))
                    mono.append(md)

            # histogram for each rank
            n = 0
            if len(mono) > 0:
                mono = np_array(mono)
                no_inference = np_array(no_inference)
                poly = np_array(poly)
                binwidth = 0.025
                bins = np_arange(0, 1.0 + binwidth, binwidth)

                mono_max_count = max(np_histogram(mono, bins=bins)[0])
                mono_weights = np_ones_like(mono) * (1.0 / mono_max_count)

                w = float(len(mono)) / (len(mono) +
                                        len(poly) + len(no_inference))
                n, b, p = ax.hist(mono, bins=bins,
                                  color=(0.0, 0.0, 1.0),
                                  alpha=0.25,
                                  weights=0.9 * w * mono_weights,
                                  bottom=i,
                                  lw=0,
                                  zorder=0)

            if len(no_inference) > 0:
                no_inference_max_count = max(
                    np_histogram(no_inference, bins=bins)[0])
                no_inference_weights = np_ones_like(
                    no_inference) * (1.0 / no_inference_max_count)

                ax.hist(no_inference, bins=bins,
                        color=(0.3, 0.3, 0.3),
                        alpha=0.25,
                        weights=0.9 * (1.0 - w) * no_inference_weights,
                        bottom=i + n,
                        lw=0,
                        zorder=0)

            if len(poly) > 0:
                poly_max_count = max(np_histogram(poly, bins=bins)[0])
                poly_weights = np_ones_like(poly) * (1.0 / poly_max_count)

                ax.hist(poly, bins=bins,
                        color=(1.0, 0.0, 0.0),
                        alpha=0.25,
                        weights=0.9 * (1.0 - w) * poly_weights,
                        bottom=i + n,
                        lw=0,
                        zorder=0)

        scatter = ax.scatter(x, y, alpha=0.5, s=48, c=c, zorder=1)

        # set plot elements
        ax.grid(color=(0.8, 0.8, 0.8), linestyle='dashed')

        ax.set_xlabel('relative distance')
        ax.set_xticks(np_arange(0, 1.05, 0.1))
        ax.set_xlim([-0.01, 1.01])

        ax.set_ylabel('rank (no. taxa)')
        ax.set_yticks(list(range(0, len(medians_for_taxa))))
        ax.set_ylim([-0.2, len(medians_for_taxa) - 0.01])
        ax.set_yticklabels(rank_labels)

        self.prettify(ax)

        # make plot interactive
        mpld3.plugins.clear(self.fig)
        mpld3.plugins.connect(
            self.fig, mpld3.plugins.PointLabelTooltip(scatter, labels=labels))
        mpld3.plugins.connect(
            self.fig, mpld3.plugins.MousePosition(fontsize=10))
        mpld3.save_html(self.fig, plot_file[0:plot_file.rfind('.')] + '.html')

        self.fig.tight_layout(pad=1)
        self.fig.savefig(plot_file, dpi=self.dpi)

    def _median_summary_outlier_file(self, phylum_rel_dists,
                                     taxa_for_dist_inference,
                                     gtdb_parent_ranks,
                                     outlier_table,
                                     rank_file,
                                     verbose_table):
        """Identify outliers relative to the median of rank distributions.
        Parameters
        ----------
        phylum_rel_dists: phylum_rel_dists[phylum][rank_index][taxon] -> relative divergences
            Relative divergence of taxon at each rank for different phylum-level rootings.
        taxa_for_dist_inference : iterable
            Taxa to considered when inferring distributions.
        gtdb_parent_ranks: d[taxon] -> string indicating parent taxa
            Parent taxa for each taxon.
        outlier_table : str
            Desired name of output table.
        rank_file : str
            Desired name of file indicating median relative distance of each rank.
        verbose_table : boolean
            Print additional columns in output table.
        """

        # determine median relative distance for each taxa
        medians_for_taxa = self.taxa_median_rd(phylum_rel_dists)

        # determine median relative distance for each rank
        median_for_rank = self.rank_median_rd(
            phylum_rel_dists, taxa_for_dist_inference)

        with open(rank_file, 'w') as fout_rank:
            median_str = []
            for rank in sorted(median_for_rank.keys()):
                median_str.append(
                    '"' + Taxonomy.rank_labels[rank] + '":' + str(median_for_rank[rank]))
            fout_rank.write('{' + ','.join(median_str) + '}\n')

        fout = open(outlier_table, 'w')
        if verbose_table:
            fout.write('Taxa\tGTDB taxonomy\tMedian distance')
            fout.write('\tMedian of rank\tMedian difference')
            fout.write('\tClosest rank\tClassifciation\n')
        else:
            fout.write(
                'Taxa\tGTDB taxonomy\tMedian distance\tMedian difference\tClosest rank\tClassification\n')

        for rank in sorted(median_for_rank.keys()):
            for clade_label, dists in medians_for_taxa[rank].items():
                dists = np_array(dists)

                taxon_median = np_median(dists)
                delta = taxon_median - median_for_rank[rank]

                closest_rank_dist = 1e10
                for test_rank, test_median in median_for_rank.items():
                    abs_dist = abs(taxon_median - test_median)
                    if abs_dist < closest_rank_dist:
                        closest_rank_dist = abs_dist
                        closest_rank = Taxonomy.rank_labels[test_rank]

                classification = "OK"
                if delta < -0.2:
                    classification = "very overclassified"
                elif delta < -0.1:
                    classification = "overclassified"
                elif delta > 0.2:
                    classification = "very underclassified"
                elif delta > 0.1:
                    classification = "underclassified"

                if verbose_table:
                    fout.write('%s\t%s\t%.2f\t%.3f\t%.3f\t%s\t%s\n' % (clade_label,
                                                                       ';'.join(
                                                                           gtdb_parent_ranks[clade_label]),
                                                                       taxon_median,
                                                                       median_for_rank[rank],
                                                                       delta,
                                                                       closest_rank,
                                                                       classification))
                else:
                    fout.write('%s\t%s\t%.3f\t%.3f\t%s\t%s\n' % (clade_label,
                                                                 ';'.join(
                                                                     gtdb_parent_ranks[clade_label]),
                                                                 taxon_median,
                                                                 delta,
                                                                 closest_rank,
                                                                 classification))
        fout.close()

    def rank_median_rd(self, phylum_rel_dists, taxa_for_dist_inference):
        """Calculate median relative divergence for each rank.

        Parameters
        ----------
        phylum_rel_dists: phylum_rel_dists[phylum][rank_index][taxon] -> relative divergences
            Relative divergence of taxon at each rank for different phylum-level rootings.
        taxa_for_dist_inference : iterable
            Taxa to considered when inferring distributions.
        """

        medians_for_taxa = self.taxa_median_rd(phylum_rel_dists)

        median_for_rank = {}
        for i, rank in enumerate(sorted(medians_for_taxa.keys())):
            v = [np_median(dists) for taxon, dists in medians_for_taxa[rank].items(
            ) if taxon in taxa_for_dist_inference]

            if v:
                median_for_rank[rank] = np_median(v)

        return median_for_rank

    def taxa_median_rd(self, phylum_rel_dists):
        """Calculate the median relative divergence for each taxon.

        Parameters
        ----------
        phylum_rel_dists: phylum_rel_dists[phylum][rank_index][taxon] -> relative divergences
            Relative divergence of taxon at each rank for different phylum-level rootings.    
        """

        medians_for_taxa = defaultdict(lambda: defaultdict(list))
        for p in phylum_rel_dists:
            for rank, d in phylum_rel_dists[p].items():
                for taxon, dist in d.items():
                    medians_for_taxa[rank][taxon].append(dist)

        return medians_for_taxa

    def _is_integer(self, s):
        """Test if a string represents an integer."""
        try:
            int(s)
            return True
        except ValueError:
            return False
