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

from gtdbtk.biolib_lite.taxonomy import Taxonomy
from gtdbtk.exceptions import GTDBTkExit


class TestTaxonomy(unittest.TestCase):

    def setUp(self):
        self.dir_tmp = tempfile.mkdtemp(prefix='gtdbtk_tmp_')

    def tearDown(self):
        shutil.rmtree(self.dir_tmp)

    def test_read(self):
        expected = {
            'GCF_005435136.1': ['d__D', 'p__P', 'c__C', 'f__F', 'g__G', 's__S1'],
            '2': ['d__D', 'p__P', 'c__C', 'f__F', 'g__G', 's__S2']
        }
        path_tax = os.path.join(self.dir_tmp, 'tax_file.tsv')
        with open(path_tax, 'w') as f:
            for k, v in expected.items():
                f.write(f'{k}\t{";".join(v)}\n')
        t = Taxonomy()
        result = t.read(path_tax)
        self.assertDictEqual(expected, result)

    def test_read_canonical(self):
        to_write = {
            'GCF_005435136.1': ['d__D', 'p__P', 'c__C', 'f__F', 'g__G', 's__S1'],
            'RS_GCF_005435135.1': ['d__D', 'p__P', 'c__C', 'f__F', 'g__G', 's__S2']
        }
        expected = {
            'G005435136': to_write['GCF_005435136.1'],
            'G005435135': to_write['RS_GCF_005435135.1'],
        }
        path_tax = os.path.join(self.dir_tmp, 'tax_file.tsv')
        with open(path_tax, 'w') as f:
            for k, v in to_write.items():
                f.write(f'{k}\t{";".join(v)}\n')
        t = Taxonomy()
        result = t.read(path_tax, canonical_ids=True)
        self.assertDictEqual(expected, result)

    def test_read_error(self):
        expected = {
            'GCF_005435136.1': ['d__D', 'p__P', 'c__C', 'f__F', 'g__G', 's__S']
        }
        path_tax = os.path.join(self.dir_tmp, 'tax_file.tsv')
        with open(path_tax, 'w') as f:
            for k, v in expected.items():
                f.write(f'{k},{";".join(v)}\n')
        t = Taxonomy()
        self.assertRaises(GTDBTkExit, t.read, path_tax)
