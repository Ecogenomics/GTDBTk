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

__prog_name__ = 'build_jplace_table.py'
__prog_desc__ = 'Precompute the per-edge RED and placement taxonomy of a reference tree, to classify placements from the jplace.'

__author__ = 'Pierre Chaumeil'
__copyright__ = 'Copyright 2026'
__credits__ = ['Pierre Chaumeil']
__license__ = 'GPL3'
__version__ = '0.0.1'
__maintainer__ = 'Pierre Chaumeil'
__email__ = 'pch@bio.aau.dk'
__status__ = 'Development'

import argparse
import os
import sys
import time

from gtdbtk.jplace_table import JplaceTableBuilder, JplaceTable, refpkg_jplace_tree


if __name__ == '__main__':
    print(__prog_name__ + ' v' + __version__ + ': ' + __prog_desc__)
    print('  by ' + __author__ + ' (' + __email__ + ')' + '\n')

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--domain', choices=['bac120', 'ar53'],
                        help='use the reference package, RED and taxonomy files and table path (mrca_red/) of GTDBTK_DATA_PATH.')
    parser.add_argument('--refpkg',
                        help='gplace reference package (e.g. pplacer/gtdb_r232_bac120.refpkg).')
    parser.add_argument('--red_file',
                        help='RED file of the same release (e.g. mrca_red/gtdbtk_r232_bac120.tsv).')
    parser.add_argument('--taxonomy_file',
                        help='GTDB taxonomy file of the same release (e.g. taxonomy/gtdb_taxonomy.tsv).')
    parser.add_argument('--output_file',
                        help='table to write (default with --domain: mrca_red/gtdbtk_<release>_<domain>.jplace_table.tsv.gz; '
                             'update the mrca_red hash in get_REF_HASHES after adding it).')
    args = parser.parse_args()

    if args.domain:
        from gtdbtk.config.common import CONFIG
        refpkg = CONFIG.GPLACE_BAC120_REF_PKG if args.domain == 'bac120' else CONFIG.GPLACE_AR53_REF_PKG
        args.refpkg = args.refpkg or os.path.join(CONFIG.GPLACE_DIR, refpkg)
        args.red_file = args.red_file or (CONFIG.MRCA_RED_BAC120 if args.domain == 'bac120' else CONFIG.MRCA_RED_AR53)
        args.taxonomy_file = args.taxonomy_file or CONFIG.TAXONOMY_FILE
        args.output_file = args.output_file or (CONFIG.JPLACE_TABLE_BAC120 if args.domain == 'bac120'
                                                else CONFIG.JPLACE_TABLE_AR53)
    if not (args.refpkg and args.red_file and args.taxonomy_file and args.output_file):
        parser.error('give --domain, or --refpkg, --red_file, --taxonomy_file and --output_file')

    try:
        sys.setrecursionlimit(1000000)
        start = time.time()
        jplace_tree = refpkg_jplace_tree(args.refpkg)
        JplaceTableBuilder(jplace_tree, args.red_file, args.taxonomy_file).write(args.output_file)
        table = JplaceTable(args.output_file)
        table.check(jplace_tree, args.red_file, args.taxonomy_file)
        print(f'Wrote {len(table.edges):,} edges of {args.refpkg} to {args.output_file} '
              f'({time.time() - start:.0f} s).')
    except SystemExit:
        print("\nControlled exit resulting from an unrecoverable error or warning.")
    except:
        print("\nUnexpected error:", sys.exc_info()[0])
        raise
