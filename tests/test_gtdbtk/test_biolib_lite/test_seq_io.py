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
import shutil
import tempfile
import unittest
from collections import defaultdict

from gtdbtk.biolib_lite.seq_io import write_fasta, read_fasta


class TestBiolibLiteCommon(unittest.TestCase):

    def setUp(self):
        self.dir_tmp = tempfile.mkdtemp(prefix='gtdbtk_tmp_')

    def tearDown(self):
        shutil.rmtree(self.dir_tmp)

    def test_write_fasta(self):
        """ Test that the fasta file is written correctly with default parameters """
        path_fasta = os.path.join(self.dir_tmp, 'fasta.fna')
        seqs = {'genome_1': 'CAGTTCAGTT', 'genome_2': 'TTAGTCA', 'genome_3': 'CAG'}
        write_fasta(seqs, path_fasta)

    def test_write_fasta_wrap(self):
        """ Test that the sequences are wrapped as specified """
        path_fasta = os.path.join(self.dir_tmp, 'fasta.fna')
        seqs = {'genome_1': 'CAGTTCAGTT', 'genome_2': 'TTAGTCA', 'genome_3': 'CAG'}
        write_fasta(seqs, path_fasta, wrap=4)

        self.assertDictEqual(seqs, read_fasta(path_fasta))

        file_content = defaultdict(list)
        with open(path_fasta, 'r') as f:
            cur_gid = None
            for line in f.readlines():
                line = line.strip()
                if line.startswith('>'):
                    cur_gid = line[1:]
                else:
                    file_content[cur_gid].append(line)

        for gid, seq_list in file_content.items():
            for cur_seq in seq_list:
                self.assertLessEqual(len(cur_seq), 4)
