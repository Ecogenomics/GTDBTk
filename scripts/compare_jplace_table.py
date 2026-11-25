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

__prog_name__ = 'compare_jplace_table.py'
__prog_desc__ = 'Compare speed and results of RED / placement taxonomy from the tog tree and from the jplace table.'

__author__ = 'Pierre Chaumeil'
__copyright__ = 'Copyright 2026'
__credits__ = ['Pierre Chaumeil']
__license__ = 'GPL3'
__version__ = '0.0.1'
__maintainer__ = 'Pierre Chaumeil'
__email__ = 'pch@bio.aau.dk'
__status__ = 'Development'

import argparse
import logging
import os
import resource
import subprocess
import sys
import time
import types

from gtdbtk.jplace_table import JplaceTable, best_placements


def peak_gb():
    """Peak RSS of this process so far (GB; ru_maxrss is in KB on Linux)."""
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6


class _ClassifyFile(object):
    """Stands in for GPlaceClassifyFile: _get_gplace_taxonomy only adds genomes to it."""

    def __init__(self):
        self.data = {}

    def add_genome(self, gid, tax):
        self.data[gid] = tax

    def write(self):
        pass


if __name__ == '__main__':
    print(__prog_name__ + ' v' + __version__ + ': ' + __prog_desc__)
    print('  by ' + __author__ + ' (' + __email__ + ')' + '\n')

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--jplace', required=True,
                        help='gplace output of a run (classify/intermediate_results/gplace/gplace.bac120.json).')
    parser.add_argument('--domain', choices=['bac120', 'ar53'], default='bac120')
    parser.add_argument('--red_file', help='RED file (default: the one of GTDBTK_DATA_PATH).')
    parser.add_argument('--taxonomy_file', help='taxonomy file (default: the one of GTDBTK_DATA_PATH).')
    parser.add_argument('--table', help='jplace table (default: the one in mrca_red/ of GTDBTK_DATA_PATH).')
    parser.add_argument('--out_dir', default='compare_jplace_table', help='tog tree and mismatch report.')
    parser.add_argument('--gplace', default='gplace', help='gplace command.')
    args = parser.parse_args()

    sys.setrecursionlimit(1000000)
    logging.basicConfig(level=logging.WARNING)
    import gtdbtk.classify as classify_mod
    from gtdbtk.classify import Classify
    from gtdbtk.tools import standardise_taxonomy

    if args.red_file:
        classify_mod.CONFIG = types.SimpleNamespace(MRCA_RED_BAC120=args.red_file, MRCA_RED_AR53=args.red_file)
    else:
        from gtdbtk.config.common import CONFIG
        args.red_file = CONFIG.MRCA_RED_BAC120 if args.domain == 'bac120' else CONFIG.MRCA_RED_AR53
    if not args.taxonomy_file:
        from gtdbtk.config.common import CONFIG
        args.taxonomy_file = CONFIG.TAXONOMY_FILE
    if not args.table:
        from gtdbtk.config.common import CONFIG
        args.table = CONFIG.JPLACE_TABLE_BAC120 if args.domain == 'bac120' else CONFIG.JPLACE_TABLE_AR53
    os.makedirs(args.out_dir, exist_ok=True)

    # ---- new path first, so its peak memory is not hidden by the tree path's
    t0 = time.time()
    jplace_tree, best = best_placements(args.jplace)
    t_read = time.time() - t0
    t0 = time.time()
    table = JplaceTable(args.table)
    table.check(jplace_tree, args.red_file, args.taxonomy_file)
    t_table = time.time() - t0
    t0 = time.time()
    new = {gid: (table.red(*place), standardise_taxonomy(table.placement_taxonomy(place[0]), args.domain))
           for gid, place in best.items()}
    t_new = time.time() - t0
    mem_new = peak_gb()

    # ---- current path: gplace tog, then Classify on the tog tree
    stub = types.SimpleNamespace(logger=logging.getLogger('timestamp'))
    tree_file = os.path.join(args.out_dir, f'{args.domain}.tog.tree')
    t0 = time.time()
    subprocess.run([args.gplace, 'tog', '-o', tree_file, args.jplace], check=True, capture_output=True)
    t_tog = time.time() - t0
    t0 = time.time()
    tree = Classify._assign_mrca_red(stub, tree_file, args.domain)
    t_red = time.time() - t0
    msa = os.path.join(args.out_dir, 'genomes.fasta')                   # _get_gplace_taxonomy reads the ids
    with open(msa, 'w') as fh:
        fh.write(''.join(f'>{gid}\nA\n' for gid in best))
    t0 = time.time()
    old_tax = Classify._get_gplace_taxonomy(stub, _ClassifyFile(), args.domain, msa, tree)
    t_tax = time.time() - t0
    mem_old = peak_gb()

    # ---- compare, genome by genome
    diffs, mismatches = [], []
    for leaf in tree.leaf_node_iter():
        gid = leaf.taxon.label
        if gid not in best:
            continue
        cur = leaf.parent_node                          # first node with a reference genome below
        while not any(x.taxon.label not in best for x in cur.leaf_iter()):
            cur = cur.parent_node
        new_red, new_tax = new[gid]
        diff = abs(new_red - cur.rel_dist)
        diffs.append(diff)
        if diff > 1e-5 or new_tax != old_tax.get(gid):
            mismatches.append((gid, *best[gid], cur.rel_dist, new_red, diff, old_tax.get(gid), new_tax))

    report = os.path.join(args.out_dir, f'{args.domain}.mismatches.tsv')
    with open(report, 'w') as fh:
        fh.write('genome\tedge_num\tdistal_length\tred_tree\tred_table\tred_diff\ttaxonomy_tree\ttaxonomy_table\n')
        for row in mismatches:
            fh.write('\t'.join(map(str, row)) + '\n')

    n = len(diffs)
    same_tax = sum(new[g][1] == old_tax.get(g) for g in best)
    print(f'{n:,} genomes ({args.domain})\n')
    print(f'tree path : gplace tog {t_tog:7.1f} s + _assign_mrca_red {t_red:7.1f} s + taxonomy {t_tax:5.1f} s '
          f'= {t_tog + t_red + t_tax:7.1f} s   (peak memory {mem_old:.2f} GB)')
    print(f'table path: read jplace {t_read:6.1f} s + read/check table {t_table:5.1f} s + RED/taxonomy {t_new:5.2f} s '
          f'= {t_read + t_table + t_new:7.1f} s   (peak memory {mem_new:.2f} GB)\n')
    print(f'RED: max |diff| {max(diffs):.2g}, mean {sum(diffs) / n:.2g}, '
          f'{sum(d > 1e-5 for d in diffs)} genome(s) > 1e-5')
    print(f'placement taxonomy identical: {same_tax:,}/{len(best):,}')
    print(f'mismatches: {len(mismatches)} -> {report}')
