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

from gtdbtk.exceptions import *
from gtdbtk.markers import Markers


class TestMarkers(unittest.TestCase):

    def setUp(self):
        self.dir_tmp = tempfile.mkdtemp(prefix='gtdbtk_tmp_')
        self.markers = Markers(cpus=10)

    def tearDown(self):
        shutil.rmtree(self.dir_tmp)
        del self.markers

    def test__apply_mask(self):
        """ Test that applying a mask to an MSA works as expected with no pruning """
        gtdb_msa = {
            'ref_1': 'DDRMIENV',
            'ref_2': 'RDRMIE-V',
        }
        user_msa = {
            'usr_1': '-DRMIENV',
        }
        min_perc_aa = 0.0
        msa_mask = os.path.join(self.dir_tmp, 'msa_mask.mask')
        with open(msa_mask, 'w') as fm:
            fm.write('11100011')

        output_seqs, pruned_seqs = self.markers._apply_mask(gtdb_msa, user_msa, msa_mask, min_perc_aa)

        expected_out = {'ref_1': 'DDRNV',
                        'ref_2': 'RDR-V',
                        'usr_1': '-DRNV'}

        self.assertDictEqual(output_seqs, expected_out)

    def test__apply_mask_prune(self):
        """ Test that applying a mask to an MSA works as expected with pruning """
        gtdb_msa = {
            'ref_1': 'DDRMIENV',
            'ref_2': 'RDRMIE-V',
        }
        user_msa = {
            'usr_1': '-DRMIENV',
            'usr_2': '---MIE-V'
        }
        min_perc_aa = 0.5
        msa_mask = os.path.join(self.dir_tmp, 'msa_mask.mask')
        with open(msa_mask, 'w') as fm:
            fm.write('11100011')

        output_seqs, pruned_seqs = self.markers._apply_mask(gtdb_msa, user_msa, msa_mask, min_perc_aa)

        expected_out = {'ref_1': 'DDRNV',
                        'ref_2': 'RDR-V',
                        'usr_1': '-DRNV'}
        expected_pruned = {'usr_2': '----V'}

        self.assertDictEqual(output_seqs, expected_out)
        self.assertDictEqual(pruned_seqs, expected_pruned)

    def test__apply_mask__mask_len_mismatch(self):
        """ Test that an exception is thrown when the mask length does not equal the MSA """
        gtdb_msa = {
            'ref_1': 'ADRMIE'
        }
        user_msa = {
            'usr_1': 'AD-MIE'
        }
        min_perc_aa = 0.1
        msa_mask = os.path.join(self.dir_tmp, 'msa_mask.mask')
        with open(msa_mask, 'w') as fm:
            fm.write('101')

        self.assertRaises(MSAMaskLengthMismatch, self.markers._apply_mask, gtdb_msa, user_msa, msa_mask, min_perc_aa)
