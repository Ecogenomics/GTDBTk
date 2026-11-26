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

"""Checks run before the ANI screen (--de_novo_species QC file, --use_jplace_table),
and cluster numbering of classify --de_novo_species. No reference data needed."""

import argparse
import logging
import os
import shutil
import tempfile
import unittest
from unittest import mock

from gtdbtk.classify import Classify
from gtdbtk.de_novo_species import GenomeClustering
from gtdbtk.exceptions import GTDBTkExit, GenomeNotInQCFile
from gtdbtk.files.classify_summary import ClassifySummaryFileRow
from gtdbtk.main import OptionsParser


class StartupTestCase(unittest.TestCase):

    def setUp(self):
        self.dir = tempfile.mkdtemp()
        self.genome_dir = os.path.join(self.dir, 'genomes')
        os.makedirs(self.genome_dir)
        for gid in ('g1', 'g2'):
            with open(os.path.join(self.genome_dir, f'{gid}.fna'), 'w') as fh:
                fh.write(f'>{gid}_1\nACGT\n')
        # no reference data needed: skip OptionsParser.__init__ (it reads the package metadata)
        self.parser = OptionsParser.__new__(OptionsParser)
        self.parser.logger = logging.getLogger('timestamp')
        self.parser.warnings = logging.getLogger('warnings')

    def tearDown(self):
        shutil.rmtree(self.dir)

    def options(self, qc_text=None, **kwargs):
        qc_path = None
        if qc_text is not None:
            qc_path = os.path.join(self.dir, 'qc.tsv')
            with open(qc_path, 'w') as fh:
                fh.write(qc_text)
        values = dict(subparser_name='classify_wf', genome_dir=self.genome_dir, batchfile=None, extension='fna',
                      out_dir=os.path.join(self.dir, 'out'), prefix='gtdbtk', cpus=1, genes=False,
                      de_novo_species=qc_path, use_jplace_table=False, gplace_max_memory=None)
        values.update(kwargs)
        return argparse.Namespace(**values)

    def parse(self, options):
        """parse_options up to the checks: the ANI screen must never be reached."""
        with mock.patch('gtdbtk.main.check_dependencies', side_effect=AssertionError('ANI screen reached')), \
                mock.patch('gtdbtk.main.get_reference_ids', return_value=[]):
            self.parser.parse_options(options)

    def test_genome_missing_from_qc_file(self):
        with self.assertRaises(GenomeNotInQCFile):
            self.parse(self.options('Name\tCompleteness\tContamination\ng1\t99\t1\n'))

    def test_bad_header(self):
        with self.assertRaisesRegex(GTDBTkExit, 'Unknown format'):
            self.parse(self.options('id\tcomp\tcont\ng1\t99\t1\ng2\t99\t1\n', subparser_name='classify'))

    def test_valid_qc_file_passes(self):
        with self.assertRaisesRegex(AssertionError, 'ANI screen reached'):
            self.parse(self.options('Name\tCompleteness\tContamination\ng1\t99\t1\ng2\t98\t1\n'))

    def test_missing_jplace_table(self):
        missing = os.path.join(self.dir, 'missing.tsv.gz')
        config_class = 'gtdbtk.config.common.CONFIG.__class__'
        with mock.patch(f'{config_class}.JPLACE_TABLE_BAC120', new_callable=mock.PropertyMock, return_value=missing), \
                mock.patch(f'{config_class}.JPLACE_TABLE_AR53', new_callable=mock.PropertyMock, return_value=missing):
            with self.assertRaisesRegex(GTDBTkExit, 'use_jplace_table'):
                self.parse(self.options(use_jplace_table=True))


class TestStandaloneNumbering(unittest.TestCase):
    """classify --de_novo_species clusters all genomes once: numbers are unique across domains."""

    def setUp(self):
        self.dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.dir)

    def test_numbers_unique_across_domains(self):
        genomes = {}
        for gid in ('ar1', 'bac1', 'bac2'):
            genomes[gid] = os.path.join(self.dir, f'{gid}.fna')
            with open(genomes[gid], 'w') as fh:
                fh.write(f'>{gid}_1\nACGT\n')
        qc_path = os.path.join(self.dir, 'qc.tsv')
        with open(qc_path, 'w') as fh:
            fh.write('Name\tCompleteness\tContamination\nar1\t99\t1\nbac1\t98\t1\nbac2\t97\t1\n')

        classify = Classify.__new__(Classify)
        classify.logger = logging.getLogger('timestamp')
        classify.gtdb_radii, classify.cpus, classify.af_threshold = None, 1, 0.5
        with mock.patch.object(GenomeClustering, 'compute_ani', return_value={}):
            (sp_cluster, _changing), path = classify._run_de_novo_species(
                qc_path, genomes, {'ar53': [], 'bac120': []}, {}, self.dir, 'test')
        # no hits: every genome is its own species, numbered 1..3 across both domains
        self.assertEqual(sorted(c['increment'] for c in sp_cluster.values()), [1, 2, 3])
        self.assertTrue(path.endswith('test.de_novo_species.tsv'))
        self.assertEqual(GenomeClustering.read_clusters(path)[0], sp_cluster)


if __name__ == '__main__':
    unittest.main()
