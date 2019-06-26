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

from gtdbtk import tools


class TestTools(unittest.TestCase):

    def test_add_ncbi_prefix(self):
        refname = 'GCF_123.1'
        self.assertEqual(tools.add_ncbi_prefix(refname), 'RS_GCF_123.1')
        refname = 'GCA_456.1'
        self.assertEqual(tools.add_ncbi_prefix(refname), 'GB_GCA_456.1')
        refname = 'genome_1'
        self.assertEqual(tools.add_ncbi_prefix(refname), refname)

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
        file = 'gtdbtk/tests/data/genomes/genome_1.fna'
        self.assertEqual(tools.sha256(
            file), '64a8a537d4d964366d7eecefc7cb806f82db5544')

    def test_file_has_checksum(self):
        """Test that the checksum read equals the file contents"""
        dir_tmp = tempfile.mkdtemp(prefix='gtdbtk_tmp_')
        try:
            path_file = os.path.join(dir_tmp, 'file.txt')
            with open(path_file, 'w') as f:
                f.write('This is a test!\nFoo\n')
            with open(path_file + '.sha256', 'w') as f:
                f.write('ac4ecc7c839524f886136fbe2e3e17ab6a806a25')
            self.assertTrue(tools.file_has_checksum(path_file, '.sha256'))
        finally:
            shutil.rmtree(dir_tmp)

    def test_file_has_checksum__mismatch(self):
        """Test that the function fails if the checksum mismatches"""
        dir_tmp = tempfile.mkdtemp(prefix='gtdbtk_tmp_')
        try:
            path_file = os.path.join(dir_tmp, 'file.txt')
            with open(path_file, 'w') as f:
                f.write('This is a test!\nFoo\n')
            with open(path_file + '.sha256', 'w') as f:
                f.write('asd123')
            self.assertFalse(tools.file_has_checksum(path_file, '.sha256'))
        finally:
            shutil.rmtree(dir_tmp)

    def test_file_has_checksum__file_error(self):
        """Test that the function fails if either file doesn't exist"""
        dir_tmp = tempfile.mkdtemp(prefix='gtdbtk_tmp_')
        try:
            path_file = os.path.join(dir_tmp, 'file.txt')
            with open(path_file, 'w') as f:
                f.write('This is a test!\nFoo\n')
            self.assertFalse(tools.file_has_checksum(path_file, '.sha256'))
            self.assertFalse(tools.file_has_checksum('/dev/null/foo', '.sha256'))
        finally:
            shutil.rmtree(dir_tmp)


if __name__ == '__main__':
    unittest.main()
