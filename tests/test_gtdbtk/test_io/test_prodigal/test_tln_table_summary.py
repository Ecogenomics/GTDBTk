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

from gtdbtk.config.output import PATH_TLN_TABLE_SUMMARY
from gtdbtk.exceptions import GTDBTkExit
from gtdbtk.io.prodigal.tln_table_summary import TlnTableSummaryFile


class TestTlnTableSummaryFile(unittest.TestCase):

    def setUp(self):
        self.dir_tmp = tempfile.mkdtemp(prefix='gtdbtk_tmp_')

    def tearDown(self):
        self.th = None
        shutil.rmtree(self.dir_tmp)

    def test___init__(self):
        tln = TlnTableSummaryFile(self.dir_tmp, 'tst')
        self.assertDictEqual(tln.genomes, {})
        self.assertEqual(os.path.join(self.dir_tmp, PATH_TLN_TABLE_SUMMARY.format(prefix='tst')), tln.path)

    def test_add_genome(self):
        tln = TlnTableSummaryFile(self.dir_tmp, 'tst')
        tln.add_genome('a', 4)
        tln.add_genome('b', 11)
        self.assertDictEqual({'a': 4, 'b': 11}, tln.genomes)

    def test_add_genome_raises_exception(self):
        tln = TlnTableSummaryFile(self.dir_tmp, 'tst')
        tln.add_genome('a', 4)
        self.assertRaises(GTDBTkExit, tln.add_genome, 'a', 11)

    def test_write(self):
        tln = TlnTableSummaryFile(self.dir_tmp, 'tst')
        tln.add_genome('a', 4)
        tln.add_genome('b', 11)
        tln.write()

        lines = set()
        with open(tln.path) as fh:
            [lines.add(x) for x in fh.readlines()]
        self.assertSetEqual({'a\t4\n', 'b\t11\n'}, lines)
