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
__version__ = '0.0.3'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'
__status__ = 'Development'

import argparse
import logging
import os
import sys
import traceback
from collections import defaultdict, Counter

import dendropy

from gtdbtk.biolib_lite.logger import colour, logger_setup
from gtdbtk.config.output import PATH_BAC120_TREE_FILE, PATH_AR53_TREE_FILE, PATH_BAC120_SUMMARY_OUT, \
    PATH_AR53_SUMMARY_OUT
from gtdbtk.exceptions import GTDBTkExit


class Translate(object):
    """Translate GTDB to NCBI classification via majority vote."""

    def __init__(self):
        """Initialization."""
        self._logger = logging.getLogger('timestamp')
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

    def run(self, gtdbtk_output_dir, ar53_metadata_file, bac120_metadata_file,
            output_file, gtdbtk_prefix):
        """Translate GTDB to NCBI classification via majority vote."""

        # Set the output directories
        if not (ar53_metadata_file or bac120_metadata_file):
            raise GTDBTkExit('You must specify at least one of --ar53_metadata_file or --bac120_metadata_file')
        ar_summary = os.path.join(gtdbtk_output_dir,
                                  PATH_AR53_SUMMARY_OUT.format(prefix=gtdbtk_prefix)) \
            if ar53_metadata_file else None
        ar_tree = os.path.join(gtdbtk_output_dir,
                               PATH_AR53_TREE_FILE.format(prefix=gtdbtk_prefix)) \
            if ar53_metadata_file else None
        bac_summary = os.path.join(gtdbtk_output_dir,
                                   PATH_BAC120_SUMMARY_OUT.format(prefix=gtdbtk_prefix)) \
            if bac120_metadata_file else None
        bac_tree = os.path.join(gtdbtk_output_dir,
                                PATH_BAC120_TREE_FILE.format(prefix=gtdbtk_prefix)) \
            if bac120_metadata_file else None

        # Create the output file directory.
        output_dir = os.path.dirname(output_file)
        if output_dir and not os.path.isdir(output_dir):
            os.makedirs(output_dir)

        # get NCBI taxonomy string for GTDB genomes and GTDB species clusters
        ncbi_taxa = {}
        ncbi_lineages = {}
        gtdb_sp_clusters = defaultdict(set)
        for domain, metadata_file in [('archaeal', ar53_metadata_file),
                                      ('bacterial', bac120_metadata_file)]:
            # Only process those domains which have been provided as an input.
            if metadata_file is None:
                continue
                
            self._logger.info(f'Processing {domain} metadata file.')
            if not os.path.exists(metadata_file):
                raise GTDBTkExit(f'File does not exist {metadata_file}')

            with open(metadata_file, 'r', encoding='utf-8') as f:
                header = f.readline().strip().split('\t')

                ncbi_taxonomy_index = header.index('ncbi_taxonomy')
                gtdb_genome_rep_index = header.index('gtdb_genome_representative')

                for line in f.readlines():
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

        self._logger.info(f'Read NCBI taxonomy for {len(ncbi_taxa):,} genomes.')
        self._logger.info(f'Identified {len(gtdb_sp_clusters):,} GTDB species clusters.')

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
                raise GTDBTkExit(f'Majority vote domain is undefined for {rep_id}')

        self._logger.info(f'Identified {len(ncbi_sp_classification):,} GTDB '
                          f'species clusters with an NCBI classification.')

        # convert GTDB classifications to NCBI classification
        with open(output_file, 'w') as fout:
            fout.write('user_genome\tGTDB classification\tNCBI classification\n')
            for domain, summary_file, tree_file in [('Archaea', ar_summary, ar_tree),
                                                    ('Bacteria', bac_summary, bac_tree)]:
                if summary_file is None or tree_file is None:
                    self._logger.warning(f'{domain} have been skipped as no metadata file was provided.')
                    continue
                if not os.path.exists(summary_file):
                    self._logger.warning(f'{domain} have been skipped as the summary file does not exist: {summary_file}')
                    continue
                if not os.path.exists(tree_file):
                    self._logger.warning(f'{domain} have been skipped as the tree file does not exist: {summary_file}')
                    continue

                self._logger.info(f'Parsing {tree_file}')
                tree = dendropy.Tree.get_from_path(tree_file,
                                                   schema='newick',
                                                   rooting='force-rooted',
                                                   preserve_underscores=True)

                # map genomes IDs to leaf nodes
                leaf_node_map = {}
                for leaf in tree.leaf_node_iter():
                    leaf_node_map[leaf.taxon.label] = leaf

                # get majority vote NCBI classification for each user genome
                self._logger.info(f'Reclassifying genomes in {summary_file}')
                with open(summary_file) as f:
                    header = f.readline().strip().split('\t')

                    gtdb_classification_index = header.index('classification')

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

        self._logger.info(f'Results have been written to: {output_file}')


if __name__ == "__main__":
    print(__prog_name__ + ' v' + __version__ + ': ' + __prog_desc__)
    print('  by ' + __author__ + ' (' + __email__ + ')' + '\n')


    def print_help():
        print(f'''
            {colour(f'...::: GTDB to NCBI Majority Vote v{__version__} :::...', ['bright'])}

  {colour('All arguments are required from:', ['underscore'])}
    {colour('--gtdbtk_output_dir', ['bright'])}
        The output directory produced by the GTDB-Tk classify workflow.
    {colour('--output_file', ['bright'])}
        The output file to write the translated taxonomy.

  {colour('At least one argument is required from:', ['underscore'])}
    {colour('--ar53_metadata_file', ['bright'])}
        The archaeal GTDB metadata file (if processing archaeal genomes).
    {colour('--bac120_metadata_file', ['bright'])}
        The bacterial GTDB metadata file (if processing bacterial genomes).

    NOTE: Metadata files are available for download from the GTDB repository
          https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/ar53_metadata.tsv
          https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/bac120_metadata.tsv

  {colour('Optional arguments:', ['underscore'])}
    {colour('--gtdbtk_prefix', ['bright'])}
        The prefix of the GTDB-Tk output files specified in --gtdbtk_output_dir.

''')


    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--gtdbtk_output_dir', required=True,
                        help='The output directory produced by the GTDB-Tk classify workflow.')
    parser.add_argument('--ar53_metadata_file', required=False, default=None,
                        help='The archaeal GTDB metadata file (if processing archaeal genomes).')
    parser.add_argument('--bac120_metadata_file', required=False, default=None,
                        help='The bacterial GTDB metadata file (if processing bacterial genomes).')
    parser.add_argument('--output_file', required=True,
                        help='The output file to write the translated taxonomy.')
    parser.add_argument('--gtdbtk_prefix', required=False, default='gtdbtk',
                        help='The prefix of the GTDB-Tk output files specified in --gtdbtk_output_dir.')

    # Parse the arguments
    if len(sys.argv) == 1:
        print_help()
        sys.exit(0)
    elif sys.argv[1] in {'-h', '--h', '-help', '--help'}:
        print_help()
        sys.exit(0)
    else:
        args = parser.parse_args()

        logger_setup(args.out_dir if hasattr(args, 'out_dir') and args.out_dir else None,
                     "gtdbtk.log", 'GTDB to NCBI majority vote', __version__, False,
                     hasattr(args, 'debug') and args.debug)
        logger = logging.getLogger('timestamp')

    try:
        p = Translate()
        p.run(args.gtdbtk_output_dir,
              args.ar53_metadata_file,
              args.bac120_metadata_file,
              args.output_file,
              args.gtdbtk_prefix)
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
        logger.error('Controlled exit resulting from an unrecoverable error or warning.')
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
