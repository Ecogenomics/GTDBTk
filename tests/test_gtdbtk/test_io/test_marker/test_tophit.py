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

from gtdbtk.config.output import PFAM_TOP_HIT_SUFFIX, TIGRFAM_TOP_HIT_SUFFIX
from gtdbtk.io.marker.tophit import Hit, TopHitPfamFile, TopHitFile, TopHitTigrFile


class TestHit(unittest.TestCase):

    def test___init__(self):
        hit = Hit('gene_id', 'hmm_id', 0.05, 100)
        self.assertEqual('gene_id', hit.gene_id)
        self.assertEqual('hmm_id', hit.hmm_id)
        self.assertAlmostEqual(0.05, hit.e_val, 2)
        self.assertAlmostEqual(100, hit.bit_score, 0)

    def test___repr__(self):
        hit = Hit('a', 'b', 0.05, 100)
        self.assertEqual(f'a b (0.05/100)', repr(hit))

    def test___eq___equal(self):
        self.assertEqual(Hit('gene_id', 'hmm_id', 1e-100, 1e-100),
                         Hit('gene_id', 'hmm_id', 1e-100, 1e-100))

    def test___eq___not_equal(self):
        self.assertNotEqual(Hit('gene_id', 'hmm_id', 0.01, 100),
                            Hit('gene+id', 'hmm_id', 0.01, 100))
        self.assertNotEqual(Hit('gene_id', 'hmm_id', 0.01, 100),
                            Hit('gene_id', 'hmm+id', 0.01, 100))
        self.assertNotEqual(Hit('gene_id', 'hmm_id', 1e-100, 100),
                            Hit('gene_id', 'hmm_id', 1e-101, 100))
        self.assertNotEqual(Hit('gene_id', 'hmm_id', 0.01, 1e-100),
                            Hit('gene_id', 'hmm_id', 0.01, 1e-101))

    def test___lt__(self):
        self.assertLess(Hit('gene_id', 'hmm_id', 1e-100, 10),
                        Hit('gene_id', 'hmm_id', 1e-100, 20))
        self.assertLess(Hit('gene_id', 'hmm_id', 1e-100, 10),
                        Hit('gene_id', 'hmm_id', 1e-200, 10))
        self.assertLess(Hit('gene_id', 'z', 1e-100, 10),
                        Hit('gene_id', 'a', 1e-100, 10))
        self.assertLess(Hit('z', 'hmm_id', 1e-100, 10),
                        Hit('a', 'hmm_id', 1e-100, 10))

    # def test___gt__(self):
    #     self.assertGreater(Hit('gene_id', 'hmm_id', 1e-100, 20),
    #                        Hit('gene_id', 'hmm_id', 1e-100, 10))
    #     self.assertGreater(Hit('gene_id', 'hmm_id', 1e-200, 10),
    #                        Hit('gene_id', 'hmm_id', 1e-100, 10))
    #     self.assertGreater(Hit('gene_id', 'a', 1e-100, 10),
    #                        Hit('gene_id', 'z', 1e-100, 10))
    #     self.assertGreater(Hit('a', 'hmm_id', 1e-100, 10),
    #                        Hit('z', 'hmm_id', 1e-100, 10))

    def test___hash__(self):
        self.assertEqual(hash(Hit('a', 'b', 0.1, 10)),
                         hash(Hit('a', 'b', 0.1, 10)))
        self.assertNotEqual(hash(Hit('a', 'b', 0.1, 10)),
                            hash(Hit('a', 'b', 0.1, 20)))
        self.assertNotEqual(hash(Hit('a', 'b', 0.1, 10)),
                            hash(Hit('a', 'b', 0.2, 10)))
        self.assertNotEqual(hash(Hit('a', 'b', 0.1, 10)),
                            hash(Hit('a', 'c', 0.1, 10)))
        self.assertNotEqual(hash(Hit('a', 'b', 0.1, 10)),
                            hash(Hit('c', 'b', 0.1, 10)))

    def test_hmm_str(self):
        hit = Hit('gene_id', 'hmm_id', 0.05, 100)
        self.assertEqual('hmm_id,0.05,100', hit.hmm_str())

    def test_ordering(self):
        expected = [Hit('a', 'a', 0.05, 50),
                    Hit('a', 'b', 0.99, 50),
                    Hit('a', 'x', 0.01, 10),
                    Hit('a', 'y', 0.01, 10)]
        test = [expected[2], expected[1], expected[0], expected[3]]
        self.assertListEqual(expected, sorted(test, reverse=True))


class TestTopHitFile(unittest.TestCase):

    def setUp(self):
        self.dir_tmp = tempfile.mkdtemp(prefix='gtdbtk_tmp_')
        self.path = os.path.join(self.dir_tmp, 'top_hit.tsv')
        self.th = TopHitFile(self.path)

    def tearDown(self):
        self.th = None
        shutil.rmtree(self.dir_tmp)

    def test_add_hit(self):
        self.th.add_hit('a', 'b', 0.1, 10)
        self.th.add_hit('a', 'b', 0.1, 20)
        self.th.add_hit('a', 'c', 0.1, 10)
        self.th.add_hit('x', 'y', 0.1, 10)
        self.th.add_hit('x', 'y', 0.2, 10)
        self.assertEqual(Hit('a', 'b', 0.1, 20), self.th.get_top_hit('a'))
        self.assertEqual(Hit('x', 'y', 0.1, 10), self.th.get_top_hit('x'))

    def test_contains_gene_id(self):
        self.th.add_hit('a', 'b', 0.1, 10)
        self.assertTrue(self.th.contains_gene_id('a'))
        self.assertFalse(self.th.contains_gene_id('b'))

    def test_contains_gene_hmm(self):
        self.th.add_hit('a', 'b', 0.1, 10)
        self.assertTrue(self.th.contains_gene_hmm('a', 'b'))
        self.assertFalse(self.th.contains_gene_hmm('a', 'c'))
        self.assertFalse(self.th.contains_gene_hmm('c', 'b'))
        self.assertFalse(self.th.contains_gene_hmm('x', 'y'))

    def test_get_top_hit(self):
        self.th.add_hit('a', 'b', 0.1, 10)
        self.th.add_hit('a', 'c', 0.1, 20)
        self.th.add_hit('a', 'd', 0.1, 5)
        self.assertEqual(Hit('a', 'c', 0.1, 20), self.th.get_top_hit('a'))
        self.assertIsNone(self.th.get_top_hit('x'))

    def test_get_hmm_hit(self):
        self.th.add_hit('a', 'b', 0.1, 10)
        self.th.add_hit('a', 'c', 0.1, 20)
        self.assertEqual(Hit('a', 'b', 0.1, 10), self.th.get_hmm_hit('a', 'b'))

    def test_write(self):
        self.th.add_hit('a', 'b', 0.1, 10)
        self.th.add_hit('x', 'y', 0.1, 10)
        self.th.add_hit('a', 'c', 0.1, 20)
        self.th.write()

        with open(self.path, 'r') as fh:
            lines = fh.readlines()
        self.assertEqual(3, len(lines))
        self.assertEqual('Gene Id\tTop hits (Family id,e-value,bitscore)\n', lines[0])
        self.assertEqual('a\tc,0.1,20;b,0.1,10\n', lines[1])
        self.assertEqual('x\ty,0.1,10\n', lines[2])

    def test_read(self):
        with open(self.path, 'w') as fh:
            fh.write('Gene Id Top hits (Family id,e-value,bitscore)\n')
            fh.write('a\tc,0.1,20;b,0.1,10\n')
            fh.write('x\ty,0.1,10\n')
        self.th.read()
        self.assertEqual(Hit('a', 'c', 0.1, 20), self.th.get_hmm_hit('a', 'c'))
        self.assertEqual(Hit('a', 'b', 0.1, 10), self.th.get_hmm_hit('a', 'b'))
        self.assertEqual(Hit('x', 'y', 0.1, 10), self.th.get_hmm_hit('x', 'y'))

    def test_iter_hits(self):
        self.th.add_hit('a', 'b', 0.1, 10)
        self.th.add_hit('x', 'y', 0.2, 20)
        expected = {('a', Hit('a', 'b', 0.1, 10)), ('x', Hit('x', 'y', 0.2, 20))}
        self.assertSetEqual(expected, set(self.th.iter_hits()))


class TestTopHitPfamFile(unittest.TestCase):

    def setUp(self):
        self.dir_tmp = tempfile.mkdtemp(prefix='gtdbtk_tmp_')
        self.path = os.path.join(self.dir_tmp, 'top_hit.tsv')
        self.th = TopHitPfamFile(self.path, 'g')

    def tearDown(self):
        self.th = None
        shutil.rmtree(self.dir_tmp)

    def test___init__(self):
        self.assertEqual(os.path.join(self.path, 'g', f'g{PFAM_TOP_HIT_SUFFIX}'), self.th.path)


class TestTopHitTigrFile(unittest.TestCase):

    def setUp(self):
        self.dir_tmp = tempfile.mkdtemp(prefix='gtdbtk_tmp_')
        self.th = TopHitTigrFile(self.dir_tmp, 'g')

    def tearDown(self):
        self.th = None
        shutil.rmtree(self.dir_tmp)

    def test___init__(self):
        self.assertEqual(os.path.join(self.dir_tmp, 'g', f'g{TIGRFAM_TOP_HIT_SUFFIX}'), self.th.path)

    def test_add_hit(self):
        self.th.add_hit('a', 'b', 0.1, 20)
        self.th.add_hit('a', 'c', 0.1, 10)
        self.th.add_hit('x', 'y', 0.1, 20)
        self.th.add_hit('x', 'z', 0.0, 20)
        expected = {('a', Hit('a', 'b', 0.1, 20)), ('x', Hit('x', 'z', 0.0, 20))}
        self.assertSetEqual(expected, set(self.th.iter_hits()))
