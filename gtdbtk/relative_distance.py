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
import sys
import logging
from collections import defaultdict, namedtuple

from biolib.taxonomy import Taxonomy
from biolib.common import is_float
from biolib.newick import parse_label

import dendropy
import mpld3

from biolib.plots.abstract_plot import AbstractPlot

from numpy import (mean as np_mean,
                   std as np_std,
                   median as np_median,
                   abs as np_abs,
                   array as np_array,
                   arange as np_arange,
                   linspace as np_linspace,
                   percentile as np_percentile,
                   ones_like as np_ones_like,
                   histogram as np_histogram)

                   
class RelativeDistance(AbstractPlot):
    """Determine relative rates of evolutionary divergence.
    
    This code is based on Phylorank: https://github.com/dparks1134/PhyloRank 
       
    """
    
    def __init__(self):
        """Initialization."""
        self.logger = logging.getLogger()
        
        Options = namedtuple('Options', 'width height tick_font_size label_font_size dpi')
        options = Options(6, 6, 12, 12, 96)

        AbstractPlot.__init__(self, options)
        
        self.dpi = 96

    def _avg_descendant_rate(self, tree):
        """Calculate average rate of divergence for each nodes in a tree.

        The average rate is the arithmetic mean of the
        branch length to all descendant taxa.

        Parameters
        ----------
        tree : Dendropy Tree
            Phylogenetic tree.

        Returns
        -------
        The following attributes are added to each node:
          mean_dist: mean distance to tips
          num_taxa: number of terminal taxa
        """

        # calculate the mean branch length to extant taxa
        for node in tree.postorder_node_iter():
            avg_div = 0
            if node.is_leaf():
                node.mean_dist = 0.0
                node.num_taxa = 1
            else:
                node.num_taxa = sum([1 for _ in node.leaf_iter()])
                for c in node.child_node_iter():
                    num_tips = c.num_taxa
                    avg_div += (float(c.num_taxa) / node.num_taxa) * (c.mean_dist + c.edge_length)

            node.mean_dist = avg_div

    def decorate_rel_dist(self, tree):
        """Calculate relative distance to each internal node.

        Parameters
        ----------
        tree : Dendropy Tree
            Phylogenetic tree.

        Returns
        -------
        The following attributes are added to each node:
          mean_dist: mean distance to tips
          num_taxa: number of terminal taxa
          rel_dists: relative distance of node between root and extant organisms
        """

        self._avg_descendant_rate(tree)

        for node in tree.preorder_node_iter():
            if node == tree.seed_node:
                node.rel_dist = 0.0
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
            rel_dists[Taxonomy.rank_index[most_specific_rank]][taxon_name] = node.rel_dist

        return rel_dists

    def read_taxa_file(taxa_file):
        """Read taxa from file.
    
        Parameters
        ----------
        taxa_file : str
            File specifying taxa to consider. One per line.
    
        Returns
        -------
        set
            Taxa.
        """
    
        taxa = set()
        for line in open(taxa_file):
            taxa.add(line.rstrip().split('\t')[0])
    
        return taxa
    
    def filter_taxa_for_dist_inference(tree, taxonomy, trusted_taxa, min_children, min_support):
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
                    valid, error_msg = Taxonomy().validate_species_name(species_name, require_full=True, require_prefix=True)
                if not valid:
                    print '[Warning] Species name %s for %s is invalid: %s' % (species_name, taxon_id, error_msg)
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
    
        # restrict taxa used for inferring distribution to those with sufficient support
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
                    # no support value, so inform user if they were trying to filter on this property
                    print '[Error] Tree does not contain support values. As such, --min_support should be set to 0.'
                    continue
    
        # restrict taxa used for inferring distribution to the trusted set
        if trusted_taxa:
            taxa_for_dist_inference = trusted_taxa.intersection(taxa_for_dist_inference)
    
        return taxa_for_dist_inference
        
    def get_phyla_lineages(tree):
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
    
####### FILE CREATION#########    
    
    def _distribution_plot(self, rel_dists, taxa_for_dist_inference, distribution_table, plot_file):
        """Create plot showing the distribution of taxa at each taxonomic rank.

        Parameters
        ----------
        rel_dists: d[rank_index][taxon] -> relative divergence
            Relative divergence of taxa at each rank.
        taxa_for_dist_inference : iterable
            Taxa to considered when inferring distributions.
        distribution_table : str
            Desired name of output table with distribution information.
        plot_file : str
            Desired name of output plot.
        """

        self.fig.clear()
        self.fig.set_size_inches(12, 6)
        ax = self.fig.add_subplot(111)
        
        
        # create normal distributions
        for i, rank in enumerate(sorted(rel_dists.keys())):
            v = [dist for taxa, dist in rel_dists[rank].iteritems() if taxa in taxa_for_dist_inference]
            if len(v) < 2:
                continue
                
            u = np_mean(v)
            rv = norm(loc=u, scale=np_std(v))
            x = np_linspace(rv.ppf(0.001), rv.ppf(0.999), 1000)
            nd = rv.pdf(x)
            # ax.plot(x, 0.75 * (nd / max(nd)) + i, 'b-', alpha=0.6, zorder=2)
            # ax.plot((u, u), (i, i + 0.5), 'b-', zorder=2)

        # create percentile and classifciation boundary lines
        percentiles = {}
        for i, rank in enumerate(sorted(rel_dists.keys())):
            v = [dist for taxa, dist in rel_dists[rank].iteritems() if taxa in taxa_for_dist_inference]
            if len(v) == 0:
                continue
                
            p10, p50, p90 = np_percentile(v, [10, 50, 90])
            ax.plot((p10, p10), (i, i + 0.25), c=(0.3, 0.3, 0.3), lw=2, zorder=2)
            ax.plot((p50, p50), (i, i + 0.5), c=(0.3, 0.3, 0.3), lw=2, zorder=2)
            ax.plot((p90, p90), (i, i + 0.25), c=(0.3, 0.3, 0.3), lw=2, zorder=2)

            for b in [-0.2, -0.1, 0.1, 0.2]:
                boundary = p50 + b
                if boundary < 1.0 and boundary > 0.0:
                    if abs(b) == 0.1:
                        c = (1.0, 0.65, 0.0)  # orange
                    else:
                        c = (1.0, 0.0, 0.0)
                    ax.plot((boundary, boundary), (i, i + 0.5), c=c, lw=2, zorder=2)

            percentiles[i] = [p10, p50, p90]

    
        # create scatter plot and results table
        fout = open(distribution_table, 'w')
        fout.write('Taxa\tRelative Distance\tP10\tMedian\tP90\tPercentile outlier\n')
        x = []
        y = []
        c = []
        labels = []
        rank_labels = []
        for i, rank in enumerate(sorted(rel_dists.keys())):
            rank_label = Taxonomy.rank_labels[rank]
            rank_labels.append(rank_label + ' (%d)' % len(rel_dists[rank]))
            
            mono = []
            poly = []
            no_inference = []
            for clade_label, dist in rel_dists[rank].iteritems():
                x.append(dist)
                y.append(i)
                labels.append(clade_label)

                if self._is_integer(clade_label.split('^')[-1]):
                    # taxa with a numerical suffix after a caret indicate 
                    # polyphyletic groups when decorated with tax2tree
                    c.append((1.0, 0.0, 0.0))
                    poly.append(dist)
                elif clade_label not in taxa_for_dist_inference:
                    c.append((0.3, 0.3, 0.3))
                    no_inference.append(dist)
                else:
                    c.append((0.0, 0.0, 1.0))
                    mono.append(dist)
            
                # report results
                v = [clade_label, dist]
                if i in percentiles:
                    p10, p50, p90 = percentiles[i]
                    percentile_outlier = not (dist >= p10 and dist <= p90)
                    v += percentiles[i] + [str(percentile_outlier)]
                else:
                    percentile_outlier = 'Insufficent data to calculate percentiles'
                    v += [-1,-1,-1] + [str(percentile_outlier)]
                
                fout.write('%s\t%.2f\t%.2f\t%.2f\t%.2f\t%s\n' % tuple(v))
        
            # histogram for each rank
            mono = np_array(mono)
            no_inference = np_array(no_inference)
            poly = np_array(poly)
            binwidth = 0.025
            bins = np_arange(0, 1.0 + binwidth, binwidth)

            d = len(mono) + len(poly) + len(no_inference)
            if d == 0:
                break
                
            w = float(len(mono)) / d
            n = 0
            if len(mono) > 0:
                mono_max_count = max(np_histogram(mono, bins=bins)[0])
                mono_weights = np_ones_like(mono) * (1.0 / mono_max_count)

                n, b, p = ax.hist(mono, bins=bins,
                          color=(0.0, 0.0, 1.0),
                          alpha=0.25,
                          weights=0.9 * w * mono_weights,
                          bottom=i,
                          lw=0,
                          zorder=0)
                      
            if len(no_inference) > 0:
                no_inference_max_count = max(np_histogram(no_inference, bins=bins)[0])
                no_inference_weights = np_ones_like(no_inference) * (1.0 / no_inference_max_count)

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
                          
        fout.close()

    
        # overlay scatter plot elements
        scatter = ax.scatter(x, y, alpha=0.5, s=48, c=c, zorder=1)

        # set plot elements
        ax.grid(color=(0.8, 0.8, 0.8), linestyle='dashed')

        ax.set_xlabel('relative distance')
        ax.set_xticks(np_arange(0, 1.05, 0.1))
        ax.set_xlim([-0.05, 1.05])

        ax.set_ylabel('rank (no. taxa)')
        ax.set_yticks(xrange(0, len(rel_dists)))
        ax.set_ylim([-0.2, len(rel_dists) - 0.01])
        ax.set_yticklabels(rank_labels)

        self.prettify(ax)

        # make plot interactive
        mpld3.plugins.clear(self.fig)
        mpld3.plugins.connect(self.fig, mpld3.plugins.PointLabelTooltip(scatter, labels=labels))
        mpld3.plugins.connect(self.fig, mpld3.plugins.MousePosition(fontsize=10))
        mpld3.save_html(self.fig, plot_file[0:plot_file.rfind('.')] + '.html')

        self.fig.tight_layout(pad=1)
        self.fig.savefig(plot_file, dpi=self.dpi)
        

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
            v = [np_median(dists) for taxon, dists in medians_for_taxa[rank].iteritems() if taxon in taxa_for_dist_inference]
            if not v:
                # not taxa at rank suitable for creating classification boundaries
                continue
            
            p10, p50, p90 = np_percentile(v, [10, 50, 90])
            ax.plot((p10, p10), (i, i + 0.25), c=(0.3, 0.3, 0.3), lw=2, zorder=2)
            ax.plot((p50, p50), (i, i + 0.5), c=(0.3, 0.3, 0.3), lw=2, zorder=2)
            ax.plot((p90, p90), (i, i + 0.25), c=(0.3, 0.3, 0.3), lw=2, zorder=2)

            for b in [-0.2, -0.1, 0.1, 0.2]:
                boundary = p50 + b
                if boundary < 1.0 and boundary > 0.0:
                    if abs(b) == 0.1:
                        c = (1.0, 0.65, 0.0)  # orange
                    else:
                        c = (1.0, 0.0, 0.0)
                    ax.plot((boundary, boundary), (i, i + 0.5), c=c, lw=2, zorder=2)

            percentiles[i] = [p10, p50, p90]

        # create scatter plot and results table
        x = []
        y = []
        c = []
        labels = []
        rank_labels = []
        for i, rank in enumerate(sorted(medians_for_taxa.keys())):
            rank_label = Taxonomy.rank_labels[rank]
            rank_labels.append(rank_label + ' (%d)' % len(medians_for_taxa[rank]))

            mono = []
            poly = []
            no_inference = []
            for clade_label, dists in medians_for_taxa[rank].iteritems():
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

                w = float(len(mono)) / (len(mono) + len(poly) + len(no_inference))
                n, b, p = ax.hist(mono, bins=bins,
                          color=(0.0, 0.0, 1.0),
                          alpha=0.25,
                          weights=0.9 * w * mono_weights,
                          bottom=i,
                          lw=0,
                          zorder=0)
                      
            if len(no_inference) > 0:
                no_inference_max_count = max(np_histogram(no_inference, bins=bins)[0])
                no_inference_weights = np_ones_like(no_inference) * (1.0 / no_inference_max_count)

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
        ax.set_yticks(xrange(0, len(medians_for_taxa)))
        ax.set_ylim([-0.2, len(medians_for_taxa) - 0.01])
        ax.set_yticklabels(rank_labels)

        self.prettify(ax)

        # make plot interactive
        mpld3.plugins.clear(self.fig)
        mpld3.plugins.connect(self.fig, mpld3.plugins.PointLabelTooltip(scatter, labels=labels))
        mpld3.plugins.connect(self.fig, mpld3.plugins.MousePosition(fontsize=10))
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
        median_for_rank = self.rank_median_rd(phylum_rel_dists, taxa_for_dist_inference)

        fout_rank = open(rank_file, 'w')
        median_str = []
        for rank in sorted(median_for_rank.keys()):
            median_str.append('"' + Taxonomy.rank_labels[rank] + '":' + str(median_for_rank[rank]))
        fout_rank.write('{' + ','.join(median_str) + '}\n')
        fout_rank.close()
            
        fout = open(outlier_table, 'w')
        if verbose_table:
            fout.write('Taxa\tGTDB taxonomy\tMedian distance')
            fout.write('\tMedian of rank\tMedian difference')
            fout.write('\tClosest rank\tClassifciation\n')
        else:
            fout.write('Taxa\tGTDB taxonomy\tMedian distance\tMedian difference\tClosest rank\tClassification\n')
        
        for rank in sorted(median_for_rank.keys()):
            for clade_label, dists in medians_for_taxa[rank].iteritems():
                dists = np_array(dists)

                taxon_median = np_median(dists)
                delta = taxon_median - median_for_rank[rank]

                closest_rank_dist = 1e10
                for test_rank, test_median in median_for_rank.iteritems():
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
                                                                       ';'.join(gtdb_parent_ranks[clade_label]),
                                                                       taxon_median,
                                                                       median_for_rank[rank],
                                                                       delta,
                                                                       closest_rank,
                                                                       classification))
                else:
                    fout.write('%s\t%s\t%.3f\t%.3f\t%s\t%s\n' % (clade_label,
                                                                   ';'.join(gtdb_parent_ranks[clade_label]),
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
                v = [np_median(dists) for taxon, dists in medians_for_taxa[rank].iteritems() if taxon in taxa_for_dist_inference]
                
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
            for rank, d in phylum_rel_dists[p].iteritems():
                for taxon, dist in d.iteritems():
                    medians_for_taxa[rank][taxon].append(dist)
                    
        return medians_for_taxa
    
    def _median_outlier_file(self, 
                                rel_dists,
                                taxa_for_dist_inference,
                                gtdb_parent_ranks, 
                                output_file):
        """Identify outliers relative to the median of rank distributions.

        Parameters
        ----------
        rel_dists: d[rank_index][taxon] -> relative divergence
            Relative divergence of taxa at each rank.
        taxa_for_dist_inference : iterable
            Taxa to considered when inferring distributions.
        gtdb_parent_ranks: d[taxon] -> string indicating parent taxa
            Parent taxa for each taxon.
        output_file : str
            Desired name of output table.
        """

        # determine median relative distance for each rank
        median_rel_dist = {}
        for rank, d in rel_dists.iteritems():
            v = [dist for taxa, dist in d.iteritems() if taxa in taxa_for_dist_inference]
            if len(v) == 0:
                continue
                
            median_rel_dist[rank] = np_median(v)

        fout = open(output_file, 'w')
        fout.write('Taxa\tGTDB taxonomy\tMedian distance\tMean difference\tClosest rank\tClassification\n')
            
        for i, rank in enumerate(sorted(rel_dists.keys())):
            for clade_label, dist in rel_dists[rank].iteritems():
                if rank in median_rel_dist:
                    delta = dist - median_rel_dist[rank]
                    closest_rank_dist = 1e10
                    for test_rank, test_median in median_rel_dist.iteritems():
                        abs_dist = abs(dist - test_median)
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

                    fout.write('%s\t%s\t%.3f\t%.3f\t%s\t%s\n' % (clade_label,
                                                                   ';'.join(gtdb_parent_ranks[clade_label]),
                                                                   dist,
                                                                   delta,
                                                                   closest_rank,
                                                                   classification))
                else:
                    fout.write('%s\t%s\t%.3f\t%.3f\t%s\t%s\n' % (clade_label,
                                                                   ';'.join(gtdb_parent_ranks[clade_label]),
                                                                   dist,
                                                                   -1,
                                                                   'NA',
                                                                   'Insufficent data to calcualte median for rank.'))
        fout.close()

    def _is_integer(self,s):
        """Test if a string represents an integer."""
        try:
            int(s)
            return True
        except ValueError:
            return False
