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

import shutil
import tempfile
import unittest

from gtdbtk.trim_msa import TrimMSA


class TestTrimMSA(unittest.TestCase):

    def setUp(self):
        self.tmp_out_dir = tempfile.mkdtemp()
        self.trim_msa = TrimMSA(cols_per_gene=0.42,
                                min_consensus=0.25,
                                max_consensus=0.95,
                                min_perc_taxa=0.50,
                                rnd_seed=42,
                                min_perc_aa=0.10,
                                out_dir=self.tmp_out_dir)

    def tearDown(self):
        shutil.rmtree(self.tmp_out_dir)

    def test_identify_valid_columns(self):
        """ Test that the expected columns are identified as valid. """
        test_seqs = {'genome_1': 'AAVWAWGG',
                     'genome_2': 'AAVW-WGA',
                     'genome_3': '--AW-GGG'}
        result = self.trim_msa.identify_valid_columns(0, 7, test_seqs)
        self.assertSetEqual(result, {2, 5})
