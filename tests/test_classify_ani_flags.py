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

"""Classify.run: the ANI screen must not run with --genes, and must not crash when
genome_domain() is skipped. No reference data or binaries needed (skani is mocked)."""

import logging
import os
import shutil
import tempfile
import unittest
from unittest import mock

from gtdbtk import classify as classify_mod
from gtdbtk.classify import Classify
from gtdbtk.config.output import PATH_FAILS


class TestClassifyANIFlags(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix='gtdbtk_tmp_')
        self.align_dir = os.path.join(self.tmp, 'align')
        self.out_dir = os.path.join(self.tmp, 'out')
        os.makedirs(self.out_dir)
        self.classify = Classify.__new__(Classify)
        self.classify.logger = logging.getLogger('timestamp')
        self.classify.cpus = 1

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def test_genes_skips_skani(self):
        """classify --genes used to log 'ANI screening will be skipped' and then run skani anyway."""
        path = os.path.join(self.align_dir, PATH_FAILS.format(prefix='gtdbtk'))
        os.makedirs(os.path.dirname(path))
        with open(path, 'w') as fh:
            fh.write('g1\tNo marker genes\n')

        with mock.patch.object(classify_mod.ANIRep, 'run_skani') as run_skani, \
                self.assertLogs('timestamp', level='WARNING') as logs:
            self.classify.run(genomes={'g1': 'g1.faa'}, align_dir=self.align_dir, out_dir=self.out_dir,
                              prefix='gtdbtk', genes=True, skip_ani_screen=False, all_failed_prodigal=True)
        run_skani.assert_not_called()
        self.assertTrue(any('--genes flag is set to True' in m for m in logs.output))

    def test_ani_screen_without_genome_domain(self):
        """all_classified_ani skips genome_domain(); the ANI screen must still find bac_ar_diff."""
        with mock.patch.object(classify_mod.ANIRep, 'run_skani', return_value={}) as run_skani:
            output_files = self.classify.run(genomes={'g1': 'g1.fna'}, align_dir=self.align_dir,
                                             out_dir=self.out_dir, prefix='gtdbtk',
                                             skip_ani_screen=False, all_classified_ani=True)
        run_skani.assert_called_once()
        self.assertEqual(output_files, {})


if __name__ == '__main__':
    unittest.main()
