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


import gtdbtk
from gtdbtk import tools


class TestTools(unittest.TestCase):

    def test_add_ncbi_prefix(self):
        refname = 'GCF_123.1'
        self.assertEqual(tools.add_ncbi_prefix(refname), 'RS_GCF_123.1')
        refname = 'GCA_456.1'
        self.assertEqual(tools.add_ncbi_prefix(refname), 'GB_GCA_456.1')

    def test_splitchunks(self):
        test_dict = {'k1': 'v1', 'k2': 'v2',
                     'k3': 'v3', 'k4': 'v4'}
        my_gen = tools.splitchunks(test_dict, 2)
        self.assertEqual(len(next(my_gen)), 2)

    def test_splitchunks_list(self):
        test_list = [1, 2, 3, 4, 5, 6]
        my_gen = tools.splitchunks_list(test_list, 2)
        self.assertEqual(len(next(my_gen)), 3)

    def test_generateTempTableName(self):
        self.assertEqual(len(tools.generateTempTableName()), 24)

    def test_merge_two_dicts(self):
        a = {'k1': 'v1', 'k2': 'v2'}
        b = {'k3': 'v3', 'k4': 'v4'}
        join_dict = tools.merge_two_dicts(a, b)
        self.assertIn('k1', join_dict)
        self.assertIn('k3', join_dict)
        self.assertEqual(len(join_dict), 4)

    def test_sha256(self):
        file = 'tests/data/genomes/genome_template.fna'
        self.assertEqual(tools.sha256(
            file), '73b7b302834824568cfcc20f06d60e57e8bc9d6a')


if __name__ == '__main__':
    unittest.main()
