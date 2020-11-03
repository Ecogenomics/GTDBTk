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


__prog_name__ = 'trim_msa.py'
__prog_desc__ = 'Randomly select a subset of columns from the MSA of each marker.'

__author__ = 'Pierre Chaumeil and Donovan Parks'
__copyright__ = 'Copyright 2017'
__credits__ = ['Pierre Chaumeil', 'Donovan Parks']
__license__ = 'GPL3'
__version__ = '0.0.4'
__maintainer__ = 'Pierre Chaumeil'
__email__ = 'p.chaumeil@uq.edu.au'
__status__ = 'Development'

import argparse
import logging
import os
import random
import sys
from collections import defaultdict, Counter

from numpy import (mean as np_mean,
                   std as np_std)

from gtdbtk.biolib_lite.logger import logger_setup
from gtdbtk.biolib_lite.seq_io import read_fasta
from gtdbtk.exceptions import MSAMarkerLengthMismatch


class TrimMSA(object):
    """Randomly select a subset of columns from the MSA of each marker."""

    def __init__(self, cols_per_gene,
                 min_perc_aa,
                 min_consensus,
                 max_consensus,
                 min_perc_taxa,
                 rnd_seed,
                 out_dir):
        """Initialization."""

        self.output_dir = out_dir
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)

        self.subset = cols_per_gene  # default: 42 * ~120 genes ~= 5,000 columns

        # only consider columns with less than this percentage of gaps
        self.max_gaps = 1.0 - min_perc_taxa

        # only consider columns where the most common amino acid is
        # between these two percent boundaries
        self.min_identical_aa = min_consensus
        self.max_identical_aa = max_consensus

        # remove genomes without sufficient number of amino acids in MSA
        self.min_perc_aa = min_perc_aa

        random.seed(rnd_seed)

        self.logger = logging.getLogger('timestamp')

    def run(self, msa_file, marker_list):
        """Randomly select a subset of columns from the MSA of each marker."""

        # read multiple sequence alignment
        self.logger.info('Reading multiple sequence alignment.')
        msa = read_fasta(msa_file, False)
        self.logger.info('Read MSA for %d genomes.' % len(msa))

        filtered_seqs, pruned_seqs = self.trim(msa, marker_list)

        self.logger.info('Removed %d taxa have amino acids in <%.1f%% of columns in filtered MSA.' % (
            len(pruned_seqs),
            self.min_perc_aa))

        # write out trimmed sequences
        with open(os.path.join(self.output_dir, "filtered_msa.faa"), 'w') as filter_file:
            for gid, seq in filtered_seqs.items():
                fasta_outstr = ">%s\n%s\n" % (gid, seq)
                filter_file.write(fasta_outstr)

        self.logger.info('Done.')

    def trim(self, msa, marker_list):
        """ Randomly select a subset of columns from the MSA of each marker.
        Args:
            msa (dict): The multiple sequence alignment to trim.
            marker_list (list): A list of markers to use for trimming (id/name/length).

        Returns:
            (dict, dict): The trimmed MSA, and any sequences which were omitted.

        Raises:
            MSAMarkerLengthMismatch: If the MSA length does not equal the length of the marker genes.
        """
        # get marker info
        self.logger.info('Reading marker info.')
        markers = []
        total_msa_len = 0
        with open(marker_list, 'r') as f:
            f.readline()
            for line in f:
                list_info = line.split("\t")
                marker_id = list_info[0]
                marker_name = '%s: %s' % (list_info[1], list_info[2])
                marker_len = int(list_info[3])
                markers.append((marker_id, marker_name, marker_len))
                total_msa_len += marker_len

        if len(list(msa.values())[0]) == total_msa_len:
            self.logger.info(f'Length of MSA and length of marker genes both equal {total_msa_len:,} columns')
        else:
            raise MSAMarkerLengthMismatch(f'Length of MSA ({len(list(msa.values())[0]):,} columns) '
                                          f'does not equal length of marker genes ({total_msa_len} columns).')

        # randomly select columns meeting filtering criteria
        self.logger.info(f'Randomly sampling {self.subset:,} columns passing filtering criteria from each marker gene.')
        mask, output_seqs = self.subsample_msa(msa, markers)

        # write mask to file
        with open(os.path.join(self.output_dir, "mask.txt"), 'w') as mask_file:
            mask_file.write(''.join([str(n) for n in mask]))

        # write subsampled MSA to file
        nbr_aa_seqs = open(os.path.join(
            self.output_dir, "genome_msa_stats.tsv"), 'w')
        nbr_aa_seqs.write(
            'Genome ID\tMSA length\tAmino acids\tAmino acids (%)\n')
        filtered_msa = {}
        pruned_seqs = {}
        for genome_id, aligned_seq in output_seqs.items():
            aa_len = sum([1 for c in aligned_seq if c.isalpha()])
            if aa_len != 0:
                aa_perc = float(aa_len) / len(aligned_seq)
            else:
                aa_perc = 0
            len_outstr = "%s\t%d\t%d\t%.2f\n" % (
                genome_id,
                len(aligned_seq),
                aa_len,
                aa_perc * 100.0)
            nbr_aa_seqs.write(len_outstr)

            if aa_perc >= self.min_perc_aa:
                filtered_msa[genome_id] = aligned_seq
            else:
                pruned_seqs[genome_id] = aligned_seq

        nbr_aa_seqs.close()

        return filtered_msa, pruned_seqs

    def identify_valid_columns(self, start, end, seqs):
        # type: (int, int, dict) -> set
        """Identify columns meeting gap and amino acid ubiquity criteria."""

        GAP_CHARS = {'-', '.', '_', '*'}
        STANDARD_AMINO_ACIDS = set('ACDEFGHIKLMNPQRSTVWY')

        gap_count = defaultdict(int)
        amino_acids = [list() for _ in range(end - start)]
        num_genomes = 0
        for seq_id, seq in seqs.items():
            num_genomes += 1
            gene_seq = seq[start:end].upper()
            for i, ch in enumerate(gene_seq):
                if ch in GAP_CHARS:
                    gap_count[i] += 1
                else:
                    amino_acids[i].append(ch)

        valid_cols = set()
        for i in range(0, end - start):
            if float(gap_count.get(i, 0)) / num_genomes <= self.max_gaps:
                c = Counter(amino_acids[i])
                if not c.most_common(1):
                    continue

                letter, count = c.most_common(1)[0]
                if letter not in STANDARD_AMINO_ACIDS:
                    self.logger.warning(
                        'Most common amino acid was not in standard alphabet: %s' % letter)

                aa_ratio = float(count) / (num_genomes - gap_count.get(i, 0))
                if self.min_identical_aa <= aa_ratio < self.max_identical_aa:
                    valid_cols.add(i)

        return valid_cols

    def subsample_msa(self, seqs, markers):
        # type: (dict, list) -> (list, dict)
        """Sample columns from each marker in multiple sequence alignment."""

        alignment_length = len(list(seqs.values())[0])
        sampled_cols = []
        start = 0
        lack_sufficient_cols = 0
        lack_cols_marker_ids = []
        avg_perc_cols = []
        for marker_id, marker_name, marker_len in markers:
            end = start + marker_len

            valid_cols = self.identify_valid_columns(start,
                                                     end,
                                                     seqs)
            assert (len(valid_cols) <= marker_len)  # sanity check

            self.logger.info('%s: S:%d, E:%d, LEN:%d, COLS:%d, PERC:%.1f' % (
                marker_name,
                start,
                end,
                marker_len,
                len(valid_cols),
                len(valid_cols) * 100.0 / marker_len))

            avg_perc_cols.append(len(valid_cols) * 100.0 / marker_len)

            if len(valid_cols) < self.subset:
                self.logger.warning(
                    'Marker has <%d columns after filtering.' % self.subset)
                lack_sufficient_cols += 1
                lack_cols_marker_ids.append(marker_id)

            offset_valid_cols = [i + start for i in valid_cols]
            sel_cols = random.sample(offset_valid_cols, min(
                self.subset, len(offset_valid_cols)))
            sampled_cols.extend(sel_cols)

            start = end

        mask = [1 if i in sampled_cols else 0 for i in range(alignment_length)]

        self.logger.info('Identified %d of %d marker genes with <%d columns for sampling:' % (
            lack_sufficient_cols,
            len(markers),
            self.subset))
        self.logger.info('%s' % ', '.join(lack_cols_marker_ids))
        self.logger.info('Marker genes had %.1f+/-%.1f%% of columns available for selection on average.' % (
            np_mean(avg_perc_cols),
            np_std(avg_perc_cols)))
        self.logger.info('Final MSA contains %d columns.' % len(sampled_cols))

        # trim columns
        output_seqs = {}
        for seq_id, seq in seqs.items():
            masked_seq = ''.join([seq[i]
                                  for i in range(0, len(mask)) if mask[i]])
            output_seqs[seq_id] = masked_seq

        return mask, output_seqs


if __name__ == '__main__':
    print(__prog_name__ + ' v' + __version__ + ': ' + __prog_desc__)
    print('  by ' + __author__ + ' (' + __email__ + ')' + '\n')

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--msa', help='unfiltered multiple sequence alignment')
    parser.add_argument(
        '--marker_list', help='file with metadata for each marker gene')
    parser.add_argument('--cols_per_gene', type=int, default=42,
                        help='maximum number of columns to retain per gene')
    parser.add_argument('--min_perc_aa', type=float, default=0.5,
                        help='filter genomes with an insufficient percentage of AA in the MSA')
    parser.add_argument('--min_consensus', type=float, default=0.25,
                        help='minimum percentage of the same amino acid required to retain column')
    parser.add_argument('--max_consensus', type=float, default=0.95,
                        help='maximum percentage of the same amino acid required to retain column')
    parser.add_argument('--min_perc_taxa', type=float, default=0.50,
                        help='minimum percentage of taxa required to retain column')
    parser.add_argument('--rnd_seed', type=int, default=None,
                        help='random seed to use for selecting columns')
    parser.add_argument('--out_dir', help='output directory')

    args = parser.parse_args()

    logger_setup(args.out_dir, "trim_msa.log", "trim_msa", __version__, False)

    try:
        p = TrimMSA(args.cols_per_gene,
                    args.min_perc_aa,
                    args.min_consensus,
                    args.max_consensus,
                    args.min_perc_taxa,
                    args.rnd_seed,
                    args.out_dir)
        p.run(args.msa, args.marker_list)
    except SystemExit:
        print("\nControlled exit resulting from an unrecoverable error or warning.")
    except:
        print("\nUnexpected error:", sys.exc_info()[0])
        raise
