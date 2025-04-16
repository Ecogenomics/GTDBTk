#!/usr/bin/env python3

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
__prog_desc__ = 'Translate GTDB to NCBI classification via majority vote.'

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2019'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__version__ = '0.2.1'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'
__status__ = 'Development'

import os
import sys
import gzip
import argparse
import logging
import traceback
from collections import defaultdict, Counter

import dendropy

from gtdbtk.biolib_lite.logger import colour, logger_setup
from gtdbtk.config.output import (PATH_BAC120_TREE_FILE,
                                  PATH_BACKBONE_BAC120_TREE_FILE,
                                  PATH_CLASS_LEVEL_BAC120_TREE_FILE,
                                  PATH_AR53_TREE_FILE,
                                  PATH_BAC120_SUMMARY_OUT,
                                  PATH_AR53_SUMMARY_OUT)
from gtdbtk.exceptions import GTDBTkExit

FAMILY_IDX = 4
SPECIES_IDX = 6


class GtdbNcbiTranslate(object):
    """Translate GTDB to NCBI classification via majority vote."""

    def __init__(self):
        """Initialization."""

        self.logger = logging.getLogger('timestamp')
        self.rank_prefix = ['d__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__']

        self.MAX_BAC_CLASS_TREES = 100

    def parse_gtdbtk_classifications(self, gtdbtk_summary_file):
        """Parse GTDB-Tk classifications."""

        gtdbtk = {}
        with open(gtdbtk_summary_file) as f:
            header = f.readline().strip().split('\t')

            gtdb_classification_idx = header.index('classification')

            for line in f:
                tokens = line.strip().split('\t')

                gid = tokens[0]
                gtdb_taxonomy = tokens[gtdb_classification_idx]
                gtdb_taxa = [t.strip()
                             for t in gtdb_taxonomy.split(';')]

                gtdbtk[gid] = gtdb_taxa

        return gtdbtk

    def get_gtdbtk_classifications(self,
                                   ar53_metadata_file,
                                   bac120_metadata_file,
                                   gtdbtk_output_dir,
                                   gtdbtk_prefix):
        """Get GTDB-Tk classification files."""

        gtdbtk_ar_assignments = {}
        if ar53_metadata_file:
            ar_summary = os.path.join(gtdbtk_output_dir,
                                      PATH_AR53_SUMMARY_OUT.format(prefix=gtdbtk_prefix))

            if os.path.exists(ar_summary):
                gtdbtk_ar_assignments = self.parse_gtdbtk_classifications(
                    ar_summary)
            else:
                self.logger.warning(
                    'Archaeal GTDB-Tk classification file does not exist.')
                self.logger.warning(
                    'Assuming there are no archaeal genomes to reclassify.')

        gtdbtk_bac_assignments = {}
        if bac120_metadata_file:
            bac_summary = os.path.join(gtdbtk_output_dir,
                                       PATH_BAC120_SUMMARY_OUT.format(prefix=gtdbtk_prefix))

            if os.path.exists(bac_summary):
                gtdbtk_bac_assignments = self.parse_gtdbtk_classifications(
                    bac_summary)
            else:
                self.logger.warning(
                    'Bacterial GTDB-Tk classification file does not exist.')
                self.logger.warning(
                    'Assuming there are no bacterial genomes to reclassify.')

        return gtdbtk_ar_assignments, gtdbtk_bac_assignments

    def get_gtdbtk_classification_trees(self,
                                        ar53_metadata_file,
                                        bac120_metadata_file,
                                        gtdbtk_output_dir,
                                        gtdbtk_prefix):
        """Get GTDB-Tk classification trees."""

        ar_sp_tree = None
        if ar53_metadata_file:
            ar_sp_tree = os.path.join(gtdbtk_output_dir,
                                      PATH_AR53_TREE_FILE.format(prefix=gtdbtk_prefix))

        bac_sp_trees = []
        bac_backbone_tree = None
        if bac120_metadata_file:
            bac_full_tree = os.path.join(gtdbtk_output_dir,
                                         PATH_BAC120_TREE_FILE.format(prefix=gtdbtk_prefix))

            if os.path.exists(bac_full_tree):
                # GTDB-Tk was run over the full tree so we only
                # need to process this single tree with one rep
                # per GTDB species cluster
                bac_sp_trees.append(bac_full_tree)
            else:
                # GTDB-Tk was run using the divide-and-conquer trees
                bac_backbone_tree = os.path.join(gtdbtk_output_dir,
                                                 PATH_BACKBONE_BAC120_TREE_FILE.format(prefix=gtdbtk_prefix))

                if os.path.exists(bac_backbone_tree):
                    for idx in range(self.MAX_BAC_CLASS_TREES):
                        bac_sp_tree = os.path.join(gtdbtk_output_dir,
                                                   PATH_CLASS_LEVEL_BAC120_TREE_FILE.format(prefix=gtdbtk_prefix, iter=idx))
                        if os.path.exists(bac_sp_tree):
                            bac_sp_trees.append(bac_sp_tree)
                else:
                    bac_backbone_tree = None

        return ar_sp_tree, bac_sp_trees, bac_backbone_tree

    def get_ncbi_descendants(self,
                             cur_node,
                             ncbi_sp_classification,
                             leaf_to_gids=None):
        """Move up tree until lineage contains at least one NCBI-defined species cluster."""

        # traverse up tree until lineage contains >=1 species with an
        # NCBI classification
        while cur_node:
            ncbi_rep_ids = set()
            for leaf in cur_node.leaf_iter():
                # leaf nodes in the backbone tree need to be expanded
                # to the set of GTDB species representatives contained
                # in a given family
                if leaf_to_gids:
                    gids = leaf_to_gids.get(leaf.taxon.label, [])
                else:
                    gids = [leaf.taxon.label]

                for gid in gids:
                    if gid in ncbi_sp_classification:
                        ncbi_rep_ids.add(gid)

            if ncbi_rep_ids:
                break

            cur_node = cur_node.parent_node

        return ncbi_rep_ids

    def parse_gtdb_metadata(self, ar53_metadata_file, bac120_metadata_file):
        """Parse GTDB metdata files to get NCBI taxonomy information and GTDB species clusters."""

        ncbi_taxa = {}
        ncbi_lineages = {}
        gtdb_sp_clusters = defaultdict(set)
        gid_to_gtdb_family = {}
        gtdb_family_to_rids = defaultdict(set)
        gtdb_sp_to_rid = {}
        for domain, metadata_file in [('archaeal', ar53_metadata_file),
                                      ('bacterial', bac120_metadata_file)]:
            # Only process those domains which have been provided as an input.
            if metadata_file is None:
                continue

            open_file = open
            if metadata_file.endswith('.gz'):
                open_file = gzip.open

            self.logger.info(f'Processing {domain} metadata file.')
            if not os.path.exists(metadata_file):
                raise GTDBTkExit(f'File does not exist: {metadata_file}')

            with open_file(metadata_file, 'rt', encoding='utf-8') as f:
                header = f.readline().strip().split('\t')

                gtdb_taxonomy_index = header.index('gtdb_taxonomy')
                ncbi_taxonomy_index = header.index('ncbi_taxonomy')
                gtdb_genome_rep_index = header.index(
                    'gtdb_genome_representative')

                for line in f.readlines():
                    tokens = line.strip().split('\t')
                    if len(tokens) <= 1:
                        # skip blank lines or ends of gzip files
                        continue

                    gid = tokens[0]
                    ncbi_taxonomy = tokens[ncbi_taxonomy_index]

                    if ncbi_taxonomy and ncbi_taxonomy != 'none':
                        ncbi_taxa[gid] = [t.strip()
                                          for t in ncbi_taxonomy.split(';')]

                        for idx, taxon in enumerate(ncbi_taxa[gid]):
                            ncbi_lineages[taxon] = ncbi_taxa[gid][0:idx + 1]
                            if idx < 6:
                                ncbi_lineages[taxon] += self.rank_prefix[idx + 1:]

                    rep_id = tokens[gtdb_genome_rep_index]
                    if gid == rep_id:
                        # genome is a GTDB representative
                        gtdb_taxonomy = tokens[gtdb_taxonomy_index]
                        gtdb_taxa = [t.strip()
                                     for t in gtdb_taxonomy.split(';')]
                        gtdb_family = gtdb_taxa[FAMILY_IDX]
                        gtdb_family_to_rids[gtdb_family].add(gid)

                        gid_to_gtdb_family[gid] = gtdb_family
                        gtdb_sp_to_rid[gtdb_taxa[SPECIES_IDX]] = rep_id

                    gtdb_sp_clusters[rep_id].add(gid)

        return ncbi_taxa, ncbi_lineages, gtdb_sp_clusters, gid_to_gtdb_family, gtdb_family_to_rids, gtdb_sp_to_rid

    def resolve_majority_vote(self, taxon_counter, num_votes):
        """Resolve majority vote taxon.

        A named taxon is considered to have a majority vote if it has >= 50%
        of the votes and not other named taxon also has 50% of the votes. Otherwise,
        there is no majority vote taxon.

        Examples:
            Case 1: >50%
                g__Bacillus = 51%, g__Metabacillus = 49% => g__Bacillus
                g__ = 51%, g__Metabacillus = 49% => no majority vote

            Case 2: no taxon at >= 50%
                g__Bacillus = 49%, g__Metabacillus = 48% => no majority vote

            Case 3: single taxon at 50%
                g__Bacillus = 50%, g__Metabacillus = 49%, ... => g__Bacillus
                g__ = 50%, g__Metabacillus = 49%, ... => no majority vote

            Case 4: two taxa at 50%
                g__Bacillus = 50%, g__Metabacillus = 50% => no majority vote
                g__Bacillus = 50%, g__ = 50% => majority vote is g__Bacillus
        """

        if num_votes == 0:
            return None

        req_counts = 0.5 * num_votes

        # case 1: top taxon has >50% of the vote
        mv_taxon, mv_count = taxon_counter.most_common(1)[0]
        if mv_count > req_counts:
            if len(mv_taxon) > 3:
                return mv_taxon
            else:
                return None

        # case 2: no taxon at >= 50%
        if mv_count < req_counts:
            return None

        # case 3: single taxon at 50%
        ((mv1_taxon, mv1_count),
         (mv2_taxon, mv2_count)) = taxon_counter.most_common(2)

        if mv1_count >= req_counts and mv2_count < req_counts:
            if len(mv_taxon) > 3:
                return mv1_taxon
            else:
                return None

        # case 4: handle case where two taxa have exactly 50%
        if mv1_count >= req_counts and mv2_count >= req_counts:
            if len(mv1_taxon) > 3 and len(mv2_taxon) <= 3:
                return mv1_taxon
            elif len(mv1_taxon) <= 3 and len(mv2_taxon) > 3:
                return mv2_taxon
            elif len(mv1_taxon) > 3 and len(mv2_taxon) > 3:
                return None

        self.logger.error('Unexpected case while resolving majority vote.')
        assert False

    def ncbi_sp_majority_vote(self, gtdb_sp_clusters, ncbi_taxa, ncbi_lineages):
        """Get NCBI majority vote classification for each GTDB species cluster."""

        ncbi_sp_classification = defaultdict(list)
        for rep_id, cluster_ids in gtdb_sp_clusters.items():
            for rank in range(6, -1, -1):
                ncbi_taxon_list = []
                for cid in cluster_ids:
                    if cid in ncbi_taxa:
                        ncbi_taxon_list.append(ncbi_taxa[cid][rank])

                mv_taxon = self.resolve_majority_vote(
                    Counter(ncbi_taxon_list),
                    len(ncbi_taxon_list))

                if mv_taxon:
                    ncbi_sp_classification[rep_id] = ncbi_lineages[mv_taxon]
                    break

            if rep_id in ncbi_sp_classification and ncbi_sp_classification[rep_id][0] == 'd__':
                raise GTDBTkExit(
                    f'Majority vote domain is undefined for {rep_id}')

        return ncbi_sp_classification

    def get_ncbi_majority_vote(self,
                               gtdb_taxa,
                               ncbi_rep_ids,
                               ncbi_sp_classification,
                               ncbi_lineages
                               ):
        """Get NCBI majority vote classification of genome."""

        # take a majority vote over species with a NCBI classification, and
        # limit taxonomic resolution to most-specific rank reported by GTDB-Tk
        ncbi_classification = []
        for rank in range(6, -1, -1):
            if len(gtdb_taxa[rank]) == 3:
                continue

            ncbi_taxon_list = []
            for rep_id in ncbi_rep_ids:
                ncbi_taxon_list.append(
                    ncbi_sp_classification[rep_id][rank])

            mv_taxon = self.resolve_majority_vote(
                Counter(ncbi_taxon_list),
                len(ncbi_taxon_list))

            if mv_taxon:
                ncbi_classification = ncbi_lineages[mv_taxon]
                break

        if not ncbi_classification:
            ncbi_classification = self.rank_prefix

        return ';'.join(ncbi_classification)

    def ncbi_majority_vote(self,
                           gtdbtk_ar_assignments,
                           ar_sp_tree,
                           gtdbtk_bac_assignments,
                           bac_sp_trees,
                           bac_backbone_tree,
                           ncbi_lineages,
                           ncbi_sp_classification,
                           gid_to_gtdb_family,
                           gtdb_family_to_rids,
                           gtdb_sp_to_rid,
                           output_file):
        """Get NCBI majority vote classification for each user genome."""

        with open(output_file, 'w') as fout:
            fout.write(
                'Genome ID\tGTDB classification\tMajority vote NCBI classification\n')

            data = [
                (gtdbtk_ar_assignments,
                 [ar_sp_tree],
                 None),

                (gtdbtk_bac_assignments,
                 bac_sp_trees,
                 bac_backbone_tree)
            ]

            for gtdbtk_assignments, sp_trees, backbone_tree in data:
                if not gtdbtk_assignments:
                    continue

                # get NCBI majority vote classification for genomes
                # placed in species-level trees
                processed_gids = set()
                for tree_file in sp_trees:
                    if not os.path.exists(tree_file):
                        # can occur since all genomes might be classified
                        # using ANI prescreening
                        continue

                    self.logger.info(f' - parsing {tree_file}')
                    tree = dendropy.Tree.get_from_path(tree_file,
                                                       schema='newick',
                                                       rooting='force-rooted',
                                                       preserve_underscores=True)

                    # map genomes IDs to leaf nodes
                    leaf_node_map = {}
                    for leaf in tree.leaf_node_iter():
                        gid = leaf.taxon.label
                        leaf_node_map[gid] = leaf

                    # get majority vote NCBI classification for each genome in tree
                    for gid, gtdb_taxa in gtdbtk_assignments.items():
                        # check if genome is in this tree
                        if gid not in leaf_node_map:
                            continue

                        processed_gids.add(gid)

                        ncbi_rep_ids = self.get_ncbi_descendants(leaf_node_map[gid],
                                                                 ncbi_sp_classification)

                        ncbi_mv = self.get_ncbi_majority_vote(
                            gtdb_taxa,
                            ncbi_rep_ids,
                            ncbi_sp_classification,
                            ncbi_lineages
                        )

                        # write out results
                        fout.write('{}\t{}\t{}\n'.format(
                            gid,
                            ';'.join(gtdb_taxa),
                            ncbi_mv))

                # get NCBI majority vote classification for
                # any genomes only placed in the backbone tree
                remaining_gids = set(gtdbtk_assignments) - processed_gids
                if len(remaining_gids) > 0 and backbone_tree:
                    self.logger.info(f' - parsing {backbone_tree}')
                    tree = dendropy.Tree.get_from_path(backbone_tree,
                                                       schema='newick',
                                                       rooting='force-rooted',
                                                       preserve_underscores=True)

                    # map genomes IDs to leaf nodes
                    leaf_node_map = {}
                    leaf_to_gids = {}
                    for leaf in tree.leaf_node_iter():
                        gid = leaf.taxon.label
                        leaf_node_map[gid] = leaf

                        # map each non-query genome in the backbone tree
                        # to the set of GTDB species representatives contained
                        # in the corresponding family
                        if gid in gid_to_gtdb_family:
                            gtdb_family = gid_to_gtdb_family[gid]
                            leaf_to_gids[gid] = gtdb_family_to_rids[gtdb_family]

                    # get majority vote NCBI classification for each genome in tree
                    for gid, gtdb_taxa in gtdbtk_assignments.items():
                        # check if genome has already been classified
                        if gid not in remaining_gids:
                            continue

                        # check if genome is in backbone tree
                        # (it will not be in any tree if classified via ANI screen)
                        if gid not in leaf_node_map:
                            continue

                        processed_gids.add(gid)

                        ncbi_rep_ids = self.get_ncbi_descendants(
                            leaf_node_map[gid],
                            ncbi_sp_classification
                        )

                        ncbi_mv = self.get_ncbi_majority_vote(
                            gtdb_taxa,
                            ncbi_rep_ids,
                            ncbi_sp_classification,
                            ncbi_lineages
                        )

                        # write out results
                        fout.write('{}\t{}\t{}\n'.format(
                            gid,
                            ';'.join(gtdb_taxa),
                            ncbi_mv))

                # get NCBI majority vote classification for genomes
                # assigned to a GTDB species cluster via ANI screening
                # (i.e. genomes not in a reference tree)
                remaining_gids = set(gtdbtk_assignments) - processed_gids
                for gid in remaining_gids:
                    gtdb_taxa = gtdbtk_assignments[gid]
                    if gtdb_taxa[0].startswith('Unclassified'):
                        ncbi_mv = gtdb_taxa
                    else:
                        gtdb_sp_rid = gtdb_sp_to_rid[gtdb_taxa[SPECIES_IDX]]
                        ncbi_mv = ncbi_sp_classification[gtdb_sp_rid]

                    fout.write('{}\t{}\t{}\n'.format(
                        gid,
                        ';'.join(gtdb_taxa),
                        ';'.join(ncbi_mv)))

    def run(self,
            gtdbtk_output_dir,
            ar53_metadata_file,
            bac120_metadata_file,
            gtdbtk_prefix,
            output_file):
        """Translate GTDB to NCBI classification via majority vote."""

        # create output file directory if required
        output_dir = os.path.dirname(output_file)
        if output_dir and not os.path.isdir(output_dir):
            os.makedirs(output_dir)

        # get GTDB-Tk classification summary files
        self.logger.info('Parsing GTDB-Tk classifications:')
        gtdbtk_ar_assignments, gtdbtk_bac_assignments = self.get_gtdbtk_classifications(
            ar53_metadata_file,
            bac120_metadata_file,
            gtdbtk_output_dir,
            gtdbtk_prefix)

        self.logger.info(
            f' - identified {len(gtdbtk_ar_assignments):,} archaeal classifications')
        self.logger.info(
            f' - identified {len(gtdbtk_bac_assignments):,} bacterial classifications')

        # get GTDB-Tk classification trees
        self.logger.info('Identifying GTDB-Tk classification trees:')
        (ar_sp_tree,
         bac_sp_trees,
         bac_backbone_tree) = self.get_gtdbtk_classification_trees(
            ar53_metadata_file,
            bac120_metadata_file,
            gtdbtk_output_dir,
            gtdbtk_prefix)

        if ar_sp_tree:
            self.logger.info(' - identified archaeal backbone tree')

        if bac_backbone_tree:
            self.logger.info(' - identified bacterial backbone tree')
        self.logger.info(
            f' - identified {len(bac_sp_trees):,} bacterial tree(s)')

        # get NCBI taxonomy information and GTDB species clusters
        self.logger.info('Parsing NCBI taxonomy from GTDB metadata files:')

        (ncbi_taxa,
         ncbi_lineages,
         gtdb_sp_clusters,
         gid_to_gtdb_family,
         gtdb_family_to_rids,
         gtdb_sp_to_rid) = self.parse_gtdb_metadata(ar53_metadata_file, bac120_metadata_file)

        self.logger.info(
            f' - read NCBI taxonomy for {len(ncbi_taxa):,} genomes')
        self.logger.info(
            f' - identified {len(gtdb_sp_clusters):,} GTDB species clusters')
        self.logger.info(
            f' - identified genomes in {len(gtdb_family_to_rids):,} GTDB families'
        )

        # get majority vote NCBI classification for each GTDB species cluster
        self.logger.info(
            'Determining NCBI majority vote classifications for GTDB species clusters.')

        ncbi_sp_classification = self.ncbi_sp_majority_vote(
            gtdb_sp_clusters,
            ncbi_taxa,
            ncbi_lineages)

        self.logger.info(
            f' - identified {len(ncbi_sp_classification):,} GTDB species clusters with an NCBI classification')

        # convert GTDB classifications to NCBI classification
        self.logger.info(
            'Determining NCBI majority vote classification for each genome:')

        self.ncbi_majority_vote(
            gtdbtk_ar_assignments,
            ar_sp_tree,
            gtdbtk_bac_assignments,
            bac_sp_trees,
            bac_backbone_tree,
            ncbi_lineages,
            ncbi_sp_classification,
            gid_to_gtdb_family,
            gtdb_family_to_rids,
            gtdb_sp_to_rid,
            output_file)

        self.logger.info(f'Results written to: {output_file}')

def main():
    print(__prog_name__ + ' v' + __version__ + ': ' + __prog_desc__)
    print('  by ' + __author__ + ' (' + __email__ + ')' + '\n')

    def print_help():
        print(f'''
                {colour(f'...::: GTDB to NCBI Majority Vote v{__version__} :::...', ['bright'])}

      {colour('Required argument:', ['underscore'])}
        {colour('--gtdbtk_output_dir', ['bright'])}
            Output directory produced by the GTDB-Tk classify workflow.
        {colour('--output_file', ['bright'])}
            Output file to write the translated taxonomy.

      {colour('At least one argument is required from:', ['underscore'])}
        {colour('--ar53_metadata_file', ['bright'])}
            GTDB archaeal metadata file (if processing archaeal genomes).
        {colour('--bac120_metadata_file', ['bright'])}
            GTDB bacterial metadata file (if processing bacterial genomes).

        NOTE: GTDB metadata files are available for download at:
              https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/ar53_metadata.tsv.gz
              https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/bac120_metadata.tsv.gz

      {colour('Optional arguments:', ['underscore'])}
        {colour('--gtdbtk_prefix', ['bright'])}
            Prefix of the GTDB-Tk output files specified in --gtdbtk_output_dir.
    ''')

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--gtdbtk_output_dir', required=True,
                        help='Output directory produced by the GTDB-Tk classify workflow.')
    parser.add_argument('--ar53_metadata_file', required=False, default=None,
                        help='The archaeal GTDB metadata file (if processing archaeal genomes).')
    parser.add_argument('--bac120_metadata_file', required=False, default=None,
                        help='The bacterial GTDB metadata file (if processing bacterial genomes).')
    parser.add_argument('--output_file', required=True,
                        help='The output file to write the translated taxonomy.')
    parser.add_argument('--gtdbtk_prefix', required=False, default='gtdbtk',
                        help='Prefix of the GTDB-Tk output files specified in --gtdbtk_output_dir.')

    # parse and sanity check arguments
    if len(sys.argv) == 1:
        print_help()
        sys.exit(0)
    elif sys.argv[1] in {'-h', '--h', '-help', '--help'}:
        print_help()
        sys.exit(0)
    else:
        logger_setup(None,
                     "gtdbtk.log",
                     'GTDB to NCBI majority vote',
                     __version__,
                     False,
                     False)
        logger = logging.getLogger('timestamp')

        args = parser.parse_args()

        if not (args.ar53_metadata_file or args.bac120_metadata_file):
            raise GTDBTkExit(
                'You must specify at least one of --ar53_metadata_file or --bac120_metadata_file')

        # check the input files exist
        for input_file in [args.ar53_metadata_file, args.bac120_metadata_file]:
            if input_file and not os.path.exists(input_file):
                logger.error(
                    f'Specified input file does not exist: {input_file}')
                sys.exit(1)

    try:
        p = GtdbNcbiTranslate()
        p.run(args.gtdbtk_output_dir,
              args.ar53_metadata_file,
              args.bac120_metadata_file,
              args.gtdbtk_prefix,
              args.output_file)
        logger.info('Done.')
    except SystemExit:
        logger.error('Controlled exit resulting from early termination.')
        sys.exit(1)
    except KeyboardInterrupt:
        logger.error('Controlled exit resulting from interrupt signal.')
        sys.exit(1)
    except GTDBTkExit as e:
        if len(str(e)) > 0:
            logger.error('{}'.format(e))
        logger.error(
            'Controlled exit resulting from an unrecoverable error or warning.')
        sys.exit(1)
    except Exception as e:
        msg = 'Uncontrolled exit resulting from an unexpected error.\n\n'
        msg += '=' * 80 + '\n'
        msg += 'EXCEPTION: {}\n'.format(type(e).__name__)
        msg += '  MESSAGE: {}\n'.format(e)
        msg += '_' * 80 + '\n\n'
        msg += traceback.format_exc()
        msg += '=' * 80
        logger.error(msg)
        sys.exit(1)


if __name__ == "__main__":
    main()