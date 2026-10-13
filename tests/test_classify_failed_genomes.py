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

"""classify when every genome failed prodigal / alignment: must finish, write the summary
and return it (it used to crash on `list > 0`). No reference data needed."""

import logging
import os
import shutil
import tempfile
import unittest

from gtdbtk.classify import Classify
from gtdbtk.config.output import PATH_FAILS, PATH_FAILED_ALIGN_GENOMES, PATH_BAC120_SUMMARY_OUT


class TestClassifyAllFailed(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix='gtdbtk_tmp_')
        self.align_dir = os.path.join(self.tmp, 'align')
        self.out_dir = os.path.join(self.tmp, 'out')
        os.makedirs(self.out_dir)
        self.prefix = 'gtdbtk'
        self.classify = Classify.__new__(Classify)
        self.classify.logger = logging.getLogger('timestamp')

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def _write(self, path_tmpl, rows):
        path = os.path.join(self.align_dir, path_tmpl.format(prefix=self.prefix))
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, 'w') as fh:
            fh.writelines(f'{gid}\t{reason}\n' for gid, reason in rows)

    def _run(self, genomes):
        return self.classify.run(genomes=genomes, align_dir=self.align_dir, out_dir=self.out_dir,
                                 prefix=self.prefix, skip_ani_screen=True, all_failed_prodigal=True)

    def test_all_failed_prodigal_does_not_crash(self):
        self._write(PATH_FAILS, [('g1', 'Prodigal failed'), ('g2', 'Prodigal failed')])
        self._write(PATH_FAILED_ALIGN_GENOMES, [('g3', 'No marker genes')])
        genomes = {'g1': 'a', 'g2': 'b', 'g3': 'c'}

        with self.assertLogs('timestamp', level='WARNING') as logs:
            output_files = self._run(genomes)

        summary = os.path.join(self.out_dir, PATH_BAC120_SUMMARY_OUT.format(prefix=self.prefix))
        self.assertEqual(output_files, {'bac120': [summary]})
        with open(summary) as fh:
            gids = [line.split('\t')[0] for line in fh.readlines()[1:]]
        self.assertEqual(sorted(set(gids)), ['g1', 'g2', 'g3'])
        self.assertIn("3 of 3 genomes have been labeled as 'Unclassified'.", '\n'.join(logs.output))

    def test_single_failed_genome_message(self):
        self._write(PATH_FAILS, [('g1', 'Prodigal failed')])
        with self.assertLogs('timestamp', level='WARNING') as logs:
            self._run({'g1': 'a'})
        self.assertIn("1 of 1 genome has been labeled as 'Unclassified'.", '\n'.join(logs.output))


if __name__ == '__main__':
    unittest.main()
