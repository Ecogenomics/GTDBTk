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

from gtdbtk.config.output import PATH_BAC120_MARKER_SUMMARY, PATH_AR122_MARKER_SUMMARY
from gtdbtk.io.marker.copy_number import CopyNumberFile, CopyNumberFileAR122, CopyNumberFileBAC120
from gtdbtk.io.marker.tophit import TopHitPfamFile, TopHitTigrFile, Hit


class TestCopyNumberFile(unittest.TestCase):

    def setUp(self):
        self.dir_tmp = tempfile.mkdtemp(prefix='gtdbtk_tmp_')

        # Build a test set of markers
        self.test_markers = {"PFAM": ["PFAM_1.hmm", "PFAM_2.hmm", "PFAM_3.hmm"],
                             "TIGRFAM": ["TIGR_1.HMM", "TIGR_2.HMM"]}

        # Build the copy number file
        self.cn_path = os.path.join(self.dir_tmp, 'cn.tsv')
        self.cn = CopyNumberFile(self.cn_path, 'test5', self.test_markers)

        # Store the hit files
        self.pfam_th = dict()
        self.tigr_th = dict()
        self.faa = dict()

        # Test case for genome_1
        """
        Single Copy: PFAM_1, TIGR_1
        Multi-Unique: PFAM_2
        Multi-Copy: None
        Missing: PFAM_3, TIGR_2
        """
        self.pfam_th['genome_1'] = TopHitPfamFile(os.path.join(self.dir_tmp), 'genome_1')
        self.pfam_th['genome_1'].add_hit('gene_a', 'PFAM_1', 0.05, 100)
        self.pfam_th['genome_1'].add_hit('gene_b', 'PFAM_2', 0.05, 200)
        self.pfam_th['genome_1'].add_hit('gene_c', 'PFAM_2', 0.05, 100)
        self.pfam_th['genome_1'].write()

        self.tigr_th['genome_1'] = TopHitTigrFile(os.path.join(self.dir_tmp), 'genome_1')
        self.tigr_th['genome_1'].add_hit('gene_x', 'TIGR_1', 0.05, 100)

        self.faa['genome_1'] = os.path.join(self.dir_tmp, 'genome_1.faa')
        with open(self.faa['genome_1'], 'w') as fh:
            fh.write('>gene_a\n')
            fh.write('VVVVVV\n')
            fh.write('>gene_b\n')
            fh.write('AAVVPP\n')
            fh.write('>gene_c\n')
            fh.write('AAVVPP\n')
            fh.write('>gene_x\n')
            fh.write('AAAAAA\n')

        # Test case for genome_2
        """
        Single Copy: PFAM_2
        Multi-Unique: TIGR_1, TIGR_2
        Multi-Copy: PFAM_1, PFAM_3
        Missing: None
        """
        self.pfam_th['genome_2'] = TopHitPfamFile(os.path.join(self.dir_tmp), 'genome_2')
        self.pfam_th['genome_2'].add_hit('gene_a', 'PFAM_2', 0.05, 100)
        self.pfam_th['genome_2'].add_hit('gene_w', 'PFAM_1', 0.05, 100)
        self.pfam_th['genome_2'].add_hit('gene_z', 'PFAM_1', 0.05, 100)
        self.pfam_th['genome_2'].add_hit('gene_y', 'PFAM_3', 0.05, 100)
        self.pfam_th['genome_2'].add_hit('gene_x', 'PFAM_3', 0.05, 100)
        self.pfam_th['genome_2'].write()

        self.tigr_th['genome_2'] = TopHitTigrFile(os.path.join(self.dir_tmp), 'genome_2')
        self.tigr_th['genome_2'].add_hit('gene_w', 'TIGR_1', 0.05, 100)
        self.tigr_th['genome_2'].add_hit('gene_x', 'TIGR_1', 0.05, 200)
        self.tigr_th['genome_2'].add_hit('gene_y', 'TIGR_2', 0.01, 100)
        self.tigr_th['genome_2'].add_hit('gene_z', 'TIGR_2', 0.05, 100)

        self.faa['genome_2'] = os.path.join(self.dir_tmp, 'genome_2.faa')
        with open(self.faa['genome_2'], 'w') as fh:
            fh.write('>gene_a\n')
            fh.write('VVVVVV\n')
            fh.write('>gene_w\n')
            fh.write('AAAVVV\n')
            fh.write('>gene_x\n')
            fh.write('AAAVVV\n')
            fh.write('>gene_y\n')
            fh.write('VVVAAA\n')
            fh.write('>gene_z\n')
            fh.write('VVVAAA\n')

        # Add the genomes
        for gid in self.faa:
            self.cn.add_genome(gid, self.faa[gid], pfam_th=self.pfam_th[gid], tigr_th=self.tigr_th[gid])

    def tearDown(self):
        self.cn = None
        self.pfam_th = None
        self.tigr_th = None
        shutil.rmtree(self.dir_tmp)

    def test_add_genome(self):
        genome_1_truth = {'unq': {'PFAM_1': {'hit': Hit('gene_a', 'PFAM_1', 0.05, 100), 'seq': 'VVVVVV'},
                                  'TIGR_1': {'hit': Hit('gene_x', 'TIGR_1', 0.05, 100), 'seq': 'AAAAAA'}
                                  },
                          'muq': {'PFAM_2': {'hit': Hit('gene_b', 'PFAM_2', 0.05, 200), 'seq': 'AAVVPP'}},
                          'mul': {},
                          'mis': {'PFAM_3': None, 'TIGR_2': None}
                          }
        genome_2_truth = {'unq': {'PFAM_2': {'hit': Hit('gene_a', 'PFAM_2', 0.05, 100), 'seq': 'VVVVVV'}},
                          'muq': {'TIGR_1': {'hit': Hit('gene_x', 'TIGR_1', 0.05, 200), 'seq': 'AAAVVV'},
                                  'TIGR_2': {'hit': Hit('gene_y', 'TIGR_2', 0.01, 100), 'seq': 'VVVAAA'}},
                          'mul': {'PFAM_1': None, 'PFAM_3': None},
                          'mis': {}}

        self.assertDictEqual(genome_1_truth, self.cn.genomes['genome_1'])
        self.assertDictEqual(genome_2_truth, self.cn.genomes['genome_2'])

    def test__extract_marker_names(self):
        true = {'PFAM_1', 'PFAM_2', 'PFAM_3', 'TIGR_1', 'TIGR_2'}
        self.assertSetEqual(true, CopyNumberFile._extract_marker_names(self.test_markers))

    def test__merge_hit_files(self):
        pfam_th = TopHitPfamFile(os.path.join(self.dir_tmp), 'genome_1')
        pfam_th.add_hit('gene_a', 'PFAM_1', 0.05, 100)
        pfam_th.add_hit('gene_b', 'PFAM_2', 0.05, 200)
        pfam_th.add_hit('gene_c', 'PFAM_2', 0.05, 100)

        tigr_th = TopHitTigrFile(os.path.join(self.dir_tmp), 'genome_1')
        tigr_th.add_hit('gene_x', 'TIGR_1', 0.05, 100)

        expected = {'TIGR_1': [Hit('gene_x', 'TIGR_1', 0.05, 100)],
                    'PFAM_1': [Hit('gene_a', 'PFAM_1', 0.05, 100)],
                    'PFAM_2': [Hit('gene_b', 'PFAM_2', 0.05, 200),
                               Hit('gene_c', 'PFAM_2', 0.05, 100)]}
        self.assertDictEqual(expected, CopyNumberFile._merge_hit_files(pfam_th, tigr_th))

    def test_get_single_copy_hits(self):
        g1_expected = {'PFAM_1': {'hit': Hit('gene_a', 'PFAM_1', 0.05, 100), 'seq': 'VVVVVV'},
                       'TIGR_1': {'hit': Hit('gene_x', 'TIGR_1', 0.05, 100), 'seq': 'AAAAAA'},
                       'PFAM_2': {'hit': Hit('gene_b', 'PFAM_2', 0.05, 200), 'seq': 'AAVVPP'}}
        g2_expected = {'PFAM_2': {'hit': Hit('gene_a', 'PFAM_2', 0.05, 100), 'seq': 'VVVVVV'},
                       'TIGR_1': {'hit': Hit('gene_x', 'TIGR_1', 0.05, 200), 'seq': 'AAAVVV'},
                       'TIGR_2': {'hit': Hit('gene_y', 'TIGR_2', 0.01, 100), 'seq': 'VVVAAA'}}

        self.assertDictEqual(g1_expected, self.cn.get_single_copy_hits('genome_1'))
        self.assertDictEqual(g2_expected, self.cn.get_single_copy_hits('genome_2'))

    def test_write(self):
        self.cn.write()

        with open(self.cn.path) as fh:
            self.assertEqual('name\tnumber_unique_genes\tnumber_multiple_genes\t'
                             'number_multiple_unique_genes\tnumber_missing_genes\t'
                             'list_unique_genes\tlist_multiple_genes\t'
                             'list_multiple_unique_genes\tlist_missing_genes\n', fh.readline())
            self.assertEqual('genome_1\t2\t0\t1\t2\t'
                             'PFAM_1,TIGR_1\t\tPFAM_2\tPFAM_3,TIGR_2\n', fh.readline())
            self.assertEqual('genome_2\t1\t2\t2\t0\t'
                             'PFAM_2\tPFAM_1,PFAM_3\tTIGR_1,TIGR_2\t\n', fh.readline())

    def test_read(self):
        self.cn.write()

        test_cn = CopyNumberFile(self.cn.path, 'test5', self.test_markers)
        test_cn.read()

        genome_1_truth = {'unq': {'PFAM_1': None, 'TIGR_1': None
                                  },
                          'muq': {'PFAM_2': None},
                          'mul': {},
                          'mis': {'PFAM_3': None, 'TIGR_2': None}
                          }
        genome_2_truth = {'unq': {'PFAM_2': None},
                          'muq': {'TIGR_1': None, 'TIGR_2': None},
                          'mul': {'PFAM_1': None, 'PFAM_3': None},
                          'mis': {}}

        self.assertDictEqual(genome_1_truth, test_cn.genomes['genome_1'])
        self.assertDictEqual(genome_2_truth, test_cn.genomes['genome_2'])


class TestCopyNumberFileAR122(unittest.TestCase):

    def setUp(self):
        self.dir_tmp = tempfile.mkdtemp(prefix='gtdbtk_tmp_')
        self.cn = CopyNumberFileAR122(self.dir_tmp, 'tst')

    def tearDown(self):
        self.cn = None
        shutil.rmtree(self.dir_tmp)

    def test___init__(self):
        path = os.path.join(self.dir_tmp, PATH_AR122_MARKER_SUMMARY.format(prefix='tst'))
        self.assertEqual(path, self.cn.path)
        self.assertEqual('ar122', self.cn.marker_set)
        self.assertDictEqual({}, self.cn.genomes)
        self.assertEqual(122, len(self.cn.marker_names))


class TestCopyNumberFileBAC120(unittest.TestCase):

    def setUp(self):
        self.dir_tmp = tempfile.mkdtemp(prefix='gtdbtk_tmp_')
        self.cn = CopyNumberFileBAC120(self.dir_tmp, 'tst')

    def tearDown(self):
        self.cn = None
        shutil.rmtree(self.dir_tmp)

    def test___init__(self):
        path = os.path.join(self.dir_tmp, PATH_BAC120_MARKER_SUMMARY.format(prefix='tst'))
        self.assertEqual(path, self.cn.path)
        self.assertEqual('bac120', self.cn.marker_set)
        self.assertDictEqual({}, self.cn.genomes)
        self.assertEqual(120, len(self.cn.marker_names))
