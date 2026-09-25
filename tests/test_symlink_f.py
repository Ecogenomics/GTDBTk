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

"""symlink_f must replace existing and dangling links when re-running in the same output directory."""

import os
import shutil
import tempfile
import unittest

from gtdbtk.tools import symlink_f


class TestSymlinkF(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix='gtdbtk_tmp_')
        self.src = os.path.join(self.tmp, 'new_target.tsv')
        with open(self.src, 'w') as fh:
            fh.write('new\n')
        self.dst = os.path.join(self.tmp, 'link.tsv')

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def _assert_points_to_src(self):
        self.assertTrue(os.path.islink(self.dst))
        self.assertEqual(os.readlink(self.dst), self.src)

    def test_new_link(self):
        symlink_f(self.src, self.dst)
        self._assert_points_to_src()

    def test_replaces_dangling_link(self):
        """Re-run after the original input was moved: the old link is broken (used to raise FileExistsError)."""
        os.symlink(os.path.join(self.tmp, 'moved_away.tsv'), self.dst)
        self.assertFalse(os.path.exists(self.dst))  # dangling
        symlink_f(self.src, self.dst)
        self._assert_points_to_src()

    def test_replaces_valid_link_and_file(self):
        other = os.path.join(self.tmp, 'old_target.tsv')
        open(other, 'w').close()
        os.symlink(other, self.dst)
        symlink_f(self.src, self.dst)
        self._assert_points_to_src()
        self.assertTrue(os.path.isfile(other))  # the old target itself is untouched

        os.remove(self.dst)
        open(self.dst, 'w').close()  # regular file
        symlink_f(self.src, self.dst)
        self._assert_points_to_src()

    def test_force_false_keeps_existing(self):
        os.symlink(os.path.join(self.tmp, 'moved_away.tsv'), self.dst)
        with self.assertRaises(FileExistsError):
            symlink_f(self.src, self.dst, force=False)

    def test_never_removes_a_directory(self):
        os.mkdir(self.dst)
        open(os.path.join(self.dst, 'keep.txt'), 'w').close()
        with self.assertRaises(FileExistsError):
            symlink_f(self.src, self.dst)
        self.assertTrue(os.path.isfile(os.path.join(self.dst, 'keep.txt')))


if __name__ == '__main__':
    unittest.main()
