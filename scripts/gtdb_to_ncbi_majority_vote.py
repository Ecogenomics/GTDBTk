#!/usr/bin/env python

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


__prog_name__ = 'gtdb_to_ncbi_majority_vote.py'
__prog_desc__ = ('Translate GTDB to NCBI classification via majority vote.')

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2019'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__version__ = '0.0.2'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'
__status__ = 'Development'

import argparse
import os
import sys
from collections import defaultdict, Counter

import dendropy


class Translate(object):
    """Translate GTDB to NCBI classification via majority vote."""

    def __init__(self):
        """Initialization."""

        self.rank_prefix = ['d__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__']

    def get_ncbi_descendants(self, user_gid, tree, leaf_node_map, ncbi_sp_classification):
        """Move up tree until lineage contains at least one NCBI-defined species cluster."""

        # traverse up tree until lineage contains >=1 species with an
        # NCBI classification
        parent = leaf_node_map[user_gid]
        while parent:
            ncbi_rep_ids = set()
            for leaf in parent.leaf_iter():
                if leaf.taxon.label in ncbi_sp_classification:
                    ncbi_rep_ids.add(leaf.taxon.label)

            if ncbi_rep_ids:
                break

            parent = parent.parent_node

        return ncbi_rep_ids

    def run(self, gtdbtk_output_dir,
            ar122_metadata_file,
            bac120_metadata_file,
            output_file,
            gtdbtk_prefix):
        """Translate GTDB to NCBI classification via majority vote."""

        # get NCBI taxonomy string for GTDB genomes and GTDB species clusters
        ncbi_taxa = {}
        ncbi_lineages = {}
        gtdb_sp_clusters = defaultdict(set)
        for metadata_file in [ar122_metadata_file, bac120_metadata_file]:
            if not os.path.exists(metadata_file):
                print('[Error] File does not exists: %s' % metadata_file)
                sys.exit(-1)

            with open(metadata_file) as f:
                header = f.readline().strip().split('\t')

                ncbi_taxonomy_index = header.index('ncbi_taxonomy')
                gtdb_genome_rep_index = header.index('gtdb_genome_representative')

                for line in f:
                    line_split = line.strip().split('\t')

                    gid = line_split[0]
                    ncbi_taxonomy = line_split[ncbi_taxonomy_index]

                    if ncbi_taxonomy and ncbi_taxonomy != 'none':
                        ncbi_taxa[gid] = [t.strip() for t in ncbi_taxonomy.split(';')]

                        for idx, taxon in enumerate(ncbi_taxa[gid]):
                            ncbi_lineages[taxon] = ncbi_taxa[gid][0:idx + 1]
                            if idx < 6:
                                ncbi_lineages[taxon] += self.rank_prefix[idx + 1:]

                    rep_id = line_split[gtdb_genome_rep_index]
                    gtdb_sp_clusters[rep_id].add(gid)

        print('Read NCBI taxonomy for %d genomes.' % len(ncbi_taxa))
        print('Identified %d GTDB species clusters.' % len(gtdb_sp_clusters))

        # get majority vote NCBI classification for each GTDB species cluster
        ncbi_sp_classification = defaultdict(list)
        for rep_id, cluster_ids in gtdb_sp_clusters.items():
            for rank in range(6, -1, -1):
                ncbi_taxon_list = []
                for cid in cluster_ids:
                    if cid in ncbi_taxa:
                        ncbi_taxon_list.append(ncbi_taxa[cid][rank])

                if len(ncbi_taxon_list) > 0:
                    counter = Counter(ncbi_taxon_list)
                    mc_taxon, mc_count = counter.most_common(1)[0]

                    if mc_count >= 0.5 * len(ncbi_taxon_list) and len(mc_taxon) > 3:
                        ncbi_sp_classification[rep_id] = ncbi_lineages[mc_taxon]
                        break

            if rep_id in ncbi_sp_classification and ncbi_sp_classification[rep_id][0] == 'd__':
                print('[Error] Majority vote domain is undefined for %s.' % rep_id)
                sys.exit(-1)

        print('Identified %d GTDB species clusters with an NCBI classification.' % len(ncbi_sp_classification))

        # convert GTDB classifications to NCBI classification
        fout = open(output_file, 'w')
        fout.write('user_genome\tGTDB classification\tNCBI classification\n')
        ar_summary = os.path.join(gtdbtk_output_dir, '%s.ar122.summary.tsv' % gtdbtk_prefix)
        ar_tree = os.path.join(gtdbtk_output_dir, '%s.ar122.classify.tree' % gtdbtk_prefix)
        bac_summary = os.path.join(gtdbtk_output_dir, '%s.bac120.summary.tsv' % gtdbtk_prefix)
        bac_tree = os.path.join(gtdbtk_output_dir, '%s.bac120.classify.tree' % gtdbtk_prefix)
        for summary_file, tree_file in [(ar_summary, ar_tree), (bac_summary, bac_tree)]:
            if not os.path.exists(summary_file):
                print('[Error] File does not exists: %s' % summary_file)
                sys.exit(-1)

            if not os.path.exists(tree_file):
                print('[Error] File does not exists: %s' % tree_file)
                sys.exit(-1)

            print('Parsing %s.' % tree_file)
            tree = dendropy.Tree.get_from_path(tree_file,
                                               schema='newick',
                                               rooting='force-rooted',
                                               preserve_underscores=True)

            # map genomes IDs to leaf nodes
            leaf_node_map = {}
            for leaf in tree.leaf_node_iter():
                leaf_node_map[leaf.taxon.label] = leaf

            # get majority vote NCBI classification for each user genome
            print('Reclassifying genomes in %s.' % summary_file)
            with open(summary_file) as f:
                header = f.readline().strip().split('\t')

                gtdb_classification_index = header.index('classification')
                fastani_ref_id_index = header.index('fastani_reference')

                for line in f:
                    line_split = line.strip().split('\t')

                    user_gid = line_split[0]
                    gtdb_taxonomy = line_split[gtdb_classification_index]
                    gtdb_taxa = [t.strip() for t in gtdb_taxonomy.split(';')]
                    gtdb_species = gtdb_taxa[6]

                    ncbi_rep_ids = self.get_ncbi_descendants(user_gid,
                                                             tree,
                                                             leaf_node_map,
                                                             ncbi_sp_classification)

                    # take a majority vote over species with a NCBI classification, and
                    # limit taxonomic resolution to most-specific rank reported by GTDB-Tk
                    ncbi_classification = []
                    for rank in range(6, -1, -1):
                        if len(gtdb_taxa[rank]) == 3:
                            continue

                        ncbi_taxon_list = []
                        for rep_id in ncbi_rep_ids:
                            ncbi_taxon_list.append(ncbi_sp_classification[rep_id][rank])

                        counter = Counter(ncbi_taxon_list)
                        mc_taxon, mc_count = counter.most_common(1)[0]

                        if mc_count >= 0.5 * len(ncbi_taxon_list) and len(mc_taxon) > 3:
                            ncbi_classification = ncbi_lineages[mc_taxon]
                            break

                    # write out results
                    fout.write('%s\t%s\t%s\n' % (
                        user_gid,
                        gtdb_taxonomy,
                        ';'.join(ncbi_classification)))

        fout.close()


if __name__ == "__main__":
    print(__prog_name__ + ' v' + __version__ + ': ' + __prog_desc__)
    print('  by ' + __author__ + ' (' + __email__ + ')' + '\n')

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--gtdbtk_output_dir', required=True,
                        help='output directory produced by GTDB-Tk classify workflow')
    parser.add_argument('--ar122_metadata_file', required=True, help='archaeal GTDB metadata file')
    parser.add_argument('--bac120_metadata_file', required=True, help='bacterial GTDB metadata file')
    parser.add_argument('--output_file', required=True, help='output file')
    parser.add_argument('--gtdbtk_prefix', required=False, default='gtdbtk', help='prefix of GTDB-Tk output files')

    args = parser.parse_args()

    try:
        p = Translate()
        p.run(args.gtdbtk_output_dir,
              args.ar122_metadata_file,
              args.bac120_metadata_file,
              args.output_file,
              args.gtdbtk_prefix)
    except SystemExit:
        print("\nControlled exit resulting from an unrecoverable error or warning.")
    except:
        print("\nUnexpected error:", sys.exc_info()[0])
        raise
