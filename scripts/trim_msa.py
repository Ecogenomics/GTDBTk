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
__prog_desc__ = 'Take randomly a subset of columns from the MSA of each marker.'

__author__ = 'Pierre Chaumeil'
__copyright__ = 'Copyright 2017'
__credits__ = ['Pierre Chaumeil']
__license__ = 'GPL3'
__version__ = '0.0.1'
__maintainer__ = 'Pierre Chaumeil'
__email__ = 'p.chaumeil@uq.edu.au'
__status__ = 'Development'

import os
import sys
import argparse
import random
from biolib.seq_io import read_fasta
from collections import defaultdict, Counter


class MSATrimmer(object):
    def __init__(self):
        """Initialization."""

        self.subset = 42
        # only consider columns with less than this percentage of gaps
        self.max_gaps = 0.10
        # only consider columns with less than this percentage of identical
        # amino acids
        self.max_indentical_aa = 0.95
        self.consensus = 0.25

    def run(self, msa, mask, marker_list, taxonomy_file, metadata_file, output):
        dict_marker = {}
        dict_genomes = read_fasta(msa, False)

        sub_list_genomes = self.selectGenomes(
            dict_genomes, taxonomy_file, metadata_file)

        with open(mask, 'r') as f:
            maskstr = f.readline()

        with open(marker_list, 'r') as f:
            f.readline()
            for line in f:
                list_info = line.split("\t")
                dict_marker[list_info[0]] = int(list_info[3])

        new_mask, output_seqs = self.trim_seqs(dict_genomes,
                                               sub_list_genomes,
                                               maskstr,
                                               dict_marker)

        if not os.path.exists(output):
            os.makedirs(output)

        # Write mask
        mask_file = open(os.path.join(output, "mask.txt"), 'w')
        mask_file.write(''.join([str(n) for n in new_mask]))
        mask_file.close()

        # Write MSA
        trimmed_file = open(os.path.join(output, "trimmed_sequences.fa"), 'w')
        nbr_aa_seqs = open(os.path.join(output, "number_AA_genomes.tsv"), 'w')
        for genome_id, aligned_seq in output_seqs.iteritems():
            fasta_outstr = ">%s\n%s\n" % (genome_id, aligned_seq)
            trimmed_file.write(fasta_outstr)
            lenaa = len(aligned_seq) - (len(aligned_seq) -
                                        len(aligned_seq.replace('-', '')))
            len_outstr = "%s\t%s\t%s\n" % (genome_id, lenaa, len(aligned_seq))
            nbr_aa_seqs.write(len_outstr)
        trimmed_file.close()
        nbr_aa_seqs.close()

    def selectGenomes(self, list_genomes, taxonomy_file, metadata_file):
        dictgenusspecies = {}
        listgid = []

        # parse metadata_file:
        dict_metadata = {}
        with open(metadata_file, 'r') as mf:
            headers = mf.readline().strip().split("\t")
            ncbi_type_strain_index = headers.index('ncbi_type_material')
            checkm_completeness_index = headers.index('checkm_completeness')
            checkm_contamination_index = headers.index('checkm_contamination')

            for line in mf:
                info = line.split('\t')
                quality = float(info[checkm_completeness_index]) - \
                    5 * float(info[checkm_contamination_index])
                tm = False
                if info[ncbi_type_strain_index] == 't':
                    tm = True

                dict_metadata[info[0]] = {"strain": tm,
                                          "quality": quality
                                          }

        with open(taxonomy_file, 'r') as f:
            for line in f:
                info = line.split("\t")
                gid = info[0]
                genusspecies = info[1].split(
                    ";")[5] + info[1].split(";")[6].strip()
                if gid in list_genomes:
                    if info[1].split(";")[6].strip() != 's__':
                        if genusspecies not in dictgenusspecies:
                            dictgenusspecies[genusspecies] = gid
                        else:
                            oldgid = dictgenusspecies.get(genusspecies)

                            if dict_metadata[oldgid]['strain'] is True and dict_metadata[gid]['strain'] is False:
                                continue
                            if dict_metadata[gid]['strain'] is True and dict_metadata[oldgid]['strain'] is False:
                                dictgenusspecies[genusspecies] = gid
                            elif dict_metadata[gid]['quality'] >= dict_metadata[oldgid]['quality']:
                                dictgenusspecies[genusspecies] = gid
                            else:
                                continue
                    else:
                        listgid.append(gid)

        for _k, v in dictgenusspecies.iteritems():
            listgid.append(v)

        return listgid

    def identify_valid_columns(self, start, end, seqs, sublistgid):
        """Identify columns meeting gap and amino acid ubiquity criteria."""

        gap_count = defaultdict(int)
        amino_acids = [list() for _ in xrange(end - start)]
        for seq_id, seq in seqs.iteritems():
            if seq_id not in sublistgid:
                continue

            gene_seq = seq[start:end]
            for i, ch in enumerate(gene_seq):
                if ch == '_' or ch == '.':
                    gap_count[i] += 1
                else:
                    amino_acids[i].append(ch)

        valid_cols = set()
        for i in xrange(0, end - start):
            if float(gap_count.get(i, 0)) / len(sublistgid) <= self.max_gaps:
                c = Counter(amino_acids[i])
                if not c.most_common(1):
                    aa_ratio = 0
                else:
                    _letter, count = c.most_common(1)[0]
                    aa_ratio = float(count) / len(sublistgid)

                if aa_ratio >= self.consensus and c.most_common(1)[0][0] != '-':
                    valid_cols.add(i)
        return valid_cols

    def trim_seqs(self, seqs, sublistgid, mask, dict_marker):
        """Trim multiple sequence alignment."""

        alignment_length = len(seqs.values()[0])
        all_coords = []
        start = 0
        for marker in sorted(dict_marker):
            end = start + dict_marker[marker]
            submask = mask[start:end]

            valid_cols = self.identify_valid_columns(
                start, end, seqs, sublistgid)

            print start, end, end - start, len(valid_cols)
            # if len(valid_cols) == 0:
            #    sys.exit()

            coords = []
            for i, presence in enumerate(submask):
                if presence == '1' and i in valid_cols:
                    coords.append(start + i)

            if len(coords) < self.subset:
                all_coords.extend(coords)
            else:
                all_coords.extend(random.sample(coords, self.subset))
            start = end

        new_mask = [
            1 if x in all_coords else 0 for x in range(alignment_length)]

        # trim columns
        output_seqs = {}
        for seq_id, seq in seqs.iteritems():
            if seq_id not in sublistgid:
                continue
            masked_seq = ''.join([seq[i]
                                  for i in xrange(0, len(mask)) if new_mask[i]])
            output_seqs[seq_id] = masked_seq

        return new_mask, output_seqs


if __name__ == '__main__':
    print __prog_name__ + ' v' + __version__ + ': ' + __prog_desc__
    print '  by ' + __author__ + ' (' + __email__ + ')' + '\n'

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--msa', help='Untrimmed multi sequence alignment.')
    parser.add_argument('--mask', help='Mask file generated by gtdb.')
    parser.add_argument(
        '--marker_list', help='File listing the metadata for each markers.')
    parser.add_argument('--taxonomy_file', help='Taxonomy file from GTDB.')
    parser.add_argument('--metadata_file', help='Metadata file from GTDB.')

    parser.add_argument('--output', help='Output directory.')

    args = parser.parse_args()

    try:
        msatrimmer = MSATrimmer()
        msatrimmer.run(args.msa, args.mask, args.marker_list,
                       args.taxonomy_file, args.metadata_file, args.output)
    except SystemExit:
        print "\nControlled exit resulting from an unrecoverable error or warning."
    except:
        print "\nUnexpected error:", sys.exc_info()[0]
        raise
