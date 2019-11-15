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

from gtdbtk.biolib_lite.common import *
from gtdbtk.biolib_lite.exceptions import BioLibFileNotFound, BioLibDirNotFound


class TestBiolibLiteCommon(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_is_float(self):
        """ Check if a string can be converted to a float """
        self.assertTrue(is_float('0.1'))
        self.assertTrue(is_float('-0.1'))
        self.assertTrue(is_float('100.0'))
        self.assertTrue(is_float('1'))
        self.assertFalse(is_float('a'))
        self.assertFalse(is_float(''))

    def test_check_file_exists(self):
        """ Check if a file exists """
        tmp_out_dir = tempfile.mkdtemp(prefix='gtdbtk_tmp_')
        try:
            invalid_path = os.path.join(tmp_out_dir, 'foo.txt')
            valid_path = os.path.join(tmp_out_dir, 'test_file.txt')
            open(valid_path, 'w').close()
            self.assertRaises(BioLibFileNotFound, check_file_exists, invalid_path)
            self.assertTrue(check_file_exists(valid_path))
        finally:
            shutil.rmtree(tmp_out_dir)

    def test_check_dir_exists(self):
        """ Check if a directory exists """
        tmp_out_dir = tempfile.mkdtemp()
        try:
            tmp_out_dir = tempfile.mkdtemp()
            self.assertTrue(check_dir_exists(tmp_out_dir))
            self.assertRaises(BioLibDirNotFound, check_dir_exists, '/tmp/path/which/doesnt/exist')
        finally:
            shutil.rmtree(tmp_out_dir)

    def test_make_sure_path_exists(self):
        """ Tests if a path is always created """
        tmp_out_dir = os.path.join(tempfile.gettempdir(), 'tmpgtdb' + next(tempfile._get_candidate_names()))
        try:
            self.assertTrue(make_sure_path_exists(''))
            self.assertTrue(make_sure_path_exists(tmp_out_dir))  # Create the directory.
            self.assertTrue(make_sure_path_exists(tmp_out_dir))  # Return True as it's already created.
            self.assertRaises(BioLibIOException, make_sure_path_exists, '/dev/null/fail')
        finally:
            shutil.rmtree(tmp_out_dir)

    def test_remove_extension(self):
        """ Test that the extension is removed correctly """
        self.assertEqual(remove_extension('foo.root.txt'), 'foo.root')
        self.assertEqual(remove_extension('foo.txt'), 'foo')
        self.assertEqual(remove_extension('foo.root.txt', 'root.txt'), 'foo')
        self.assertEqual(remove_extension('foo'), 'foo')
