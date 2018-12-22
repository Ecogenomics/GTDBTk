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
__version__ = '0.0.3'
__maintainer__ = 'Pierre Chaumeil'
__email__ = 'p.chaumeil@uq.edu.au'
__status__ = 'Development'

import os
import sys
import argparse
import random
import logging
from biolib.seq_io import read_fasta
from biolib.logger import logger_setup
from collections import defaultdict, Counter

from numpy import (mean as np_mean,
                    std as np_std)


class MSATrimmer(object):
    """Randomly select a subset of columns from the MSA of each marker."""
    
    def __init__(self, output_dir):
        """Initialization."""
        
        self.output_dir = output_dir
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)

        self.subset = 42    # 42 * ~120 genes ~= 5,000 columns
                            # also, meaning of life so seems good
        
        # only consider columns with less than this percentage of gaps
        self.max_gaps = 0.50
        
        # only consider columns where the most common amino acid is
        # between these two percent boundaries
        self.max_identical_aa = 0.95
        self.min_identical_aa = 0.25
        
        logger_setup(output_dir, "trim_msa.log", "trim_msa", __version__, False)
        self.logger = logging.getLogger('timestamp')

    def run(self, msa_file, marker_list):
        """Randomly select a subset of columns from the MSA of each marker."""
        
        # read multiple sequence alignment
        self.logger.info('Reading multiple sequence alignment.')
        msa = read_fasta(msa_file, False)
        self.logger.info('Read MSA for %d genomes.' % len(msa))
        
        # get marker info
        self.logger.info('Reading marker info.')
        markers = {}
        total_msa_len = 0
        with open(marker_list, 'r') as f:
            f.readline()
            for line in f:
                list_info = line.split("\t")
                marker_id = list_info[0]
                marker_name = '%s: %s' % (list_info[1], list_info[2])
                marker_len = int(list_info[3])
                markers[marker_id] = (marker_name, marker_len)
                total_msa_len += marker_len
                
        if len(msa.values()[0]) == total_msa_len:
            self.logger.info('Length of MSA and length of marker genes both equal %d columns' % total_msa_len)
        else:
            self.logger.error('Length of MSA (%d columns) does not equal length of marker genes (%d columns).' % (
                                    len(msa.values()[0]), 
                                    total_msa_len))
            sys.exit(-1)

        # randomly select columns meeting filtering criteria
        self.logger.info('Randomly sampling %d columns passing filtering criteria from each marker gene.' % self.subset)
        mask, output_seqs = self.subsample_msa(msa, markers)

        # write mask to file
        mask_file = open(os.path.join(self.output_dir, "mask.txt"), 'w')
        mask_file.write(''.join([str(n) for n in mask]))
        mask_file.close()

        # write subsampled MSA to file
        trimmed_file = open(os.path.join(self.output_dir, "trimmed_sequences.faa"), 'w')
        nbr_aa_seqs = open(os.path.join(self.output_dir, "genome_msa_stats.tsv"), 'w')
        nbr_aa_seqs.write('Genome ID\tMSA length\tAmino acids\tAmino acids (%)\n')
        for genome_id, aligned_seq in output_seqs.iteritems():
            fasta_outstr = ">%s\n%s\n" % (genome_id, aligned_seq)
            trimmed_file.write(fasta_outstr)
            lenaa = len(aligned_seq) - (len(aligned_seq) -
                                        len(aligned_seq.replace('-', '')))
            len_outstr = "%s\t%d\t%d\t%.2f\n" % (
                            genome_id, 
                            len(aligned_seq), 
                            lenaa, 
                            lenaa*100.0/len(aligned_seq))
            nbr_aa_seqs.write(len_outstr)
        trimmed_file.close()
        nbr_aa_seqs.close()
        
        self.logger.info('Done.')

    def identify_valid_columns(self, start, end, seqs):
        """Identify columns meeting gap and amino acid ubiquity criteria."""
        
        GAP_CHARS = set(['-', '.', '_', '*'])
        STANDARD_AMINO_ACIDS = set('ACDEFGHIKLMNPQRSTVWY')
        
        gap_count = defaultdict(int)
        amino_acids = [list() for _ in xrange(end - start)]
        num_genomes = 0
        for seq_id, seq in seqs.iteritems():
            num_genomes += 1
            gene_seq = seq[start:end].upper()
            for i, ch in enumerate(gene_seq):
                if ch in GAP_CHARS:
                    gap_count[i] += 1
                else:
                    amino_acids[i].append(ch)

        valid_cols = set()
        for i in xrange(0, end - start):
            if float(gap_count.get(i, 0)) / num_genomes <= self.max_gaps:
                c = Counter(amino_acids[i])
                if not c.most_common(1):
                    continue
                
                letter, count = c.most_common(1)[0]
                if letter not in STANDARD_AMINO_ACIDS:
                    self.logger.error('Most common amino acid was not in standard alphabet: %s' % letter)
                    sys.exit(-1)

                aa_ratio = float(count) / (num_genomes - gap_count.get(i, 0))
                if self.min_identical_aa <= aa_ratio <= self.max_identical_aa:
                    valid_cols.add(i)

        return valid_cols

    def subsample_msa(self, seqs, markers):
        """Sample columns from each marker in multiple sequence alignment."""

        alignment_length = len(seqs.values()[0])
        sampled_cols = []
        start = 0
        lack_sufficient_cols = 0
        lack_cols_marker_ids = []
        avg_perc_cols = []
        for marker_id in markers:
            marker_name, marker_len = markers[marker_id]
            end = start + marker_len

            valid_cols = self.identify_valid_columns(start, 
                                                        end, 
                                                        seqs)
            assert(len(valid_cols) <= marker_len) # sanity check
            
            self.logger.info('%s: S:%d, E:%d, LEN:%d, COLS:%d, PERC:%.1f' % (
                                marker_name, 
                                start, 
                                end, 
                                marker_len, 
                                len(valid_cols),
                                len(valid_cols)*100.0/marker_len))
                                
            avg_perc_cols.append(len(valid_cols)*100.0/marker_len)
                                
            if len(valid_cols) < self.subset:
                self.logger.warning('Marker has <%d columns after filtering.' % self.subset)
                lack_sufficient_cols += 1
                lack_cols_marker_ids.append(marker_id)
                
            offset_valid_cols = [i+start for i in valid_cols]
            sampled_cols.extend(random.sample(offset_valid_cols, min(self.subset, len(offset_valid_cols))))
            
            start = end
 
        mask = [1 if i in sampled_cols else 0 for i in range(alignment_length)]
        
        self.logger.info('Identified %d of %d marker genes with <%d columns for sampling:' % (
                            lack_sufficient_cols, 
                            len(markers),
                            self.subset))
        self.logger.info('  %s' % ', '.join(lack_cols_marker_ids))
        self.logger.info('Marker genes had %.1f+/-%.1f%% of columns available for selection on average.' % (
                            np_mean(avg_perc_cols),
                            np_std(avg_perc_cols)))
        self.logger.info('Final MSA contains %d columns.' % len(sampled_cols))

        # trim columns
        output_seqs = {}
        for seq_id, seq in seqs.iteritems():
            masked_seq = ''.join([seq[i]
                                  for i in xrange(0, len(mask)) if mask[i]])
            output_seqs[seq_id] = masked_seq

        return mask, output_seqs

if __name__ == '__main__':
    print __prog_name__ + ' v' + __version__ + ': ' + __prog_desc__
    print '  by ' + __author__ + ' (' + __email__ + ')' + '\n'

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--msa', help='unfiltered multiple sequence alignment')
    parser.add_argument('--marker_list', help='file with metadata for each marker gene')
    parser.add_argument('--output', help='output directory')

    args = parser.parse_args()

    try:
        p = MSATrimmer(args.output)
        p.run(args.msa, args.marker_list)
    except SystemExit:
        print "\nControlled exit resulting from an unrecoverable error or warning."
    except:
        print "\nUnexpected error:", sys.exc_info()[0]
        raise
