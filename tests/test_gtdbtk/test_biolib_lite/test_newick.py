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

import unittest

from gtdbtk.biolib_lite.newick import *


class TestBioLibNewick(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_parse_label(self):
        """ Test that the information obtained from internal newick nodes is consistent """
        self.assertTupleEqual(parse_label('1.0:o__something|foo'), (1.0, 'o__something', 'foo'))
        self.assertTupleEqual(parse_label('1.0:o__something'), (1.0, 'o__something', None))
        self.assertTupleEqual(parse_label('1.0'), (1.0, None, None))
