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

"""Tests for the --de_novo_species rows in the classify summary file and for
skani triangle parsing. No reference data or third-party binaries are needed."""

import logging
import shutil
import tempfile
import unittest
from unittest import mock

from gtdbtk.classify import Classify
from gtdbtk.external.skani import SkANI
from gtdbtk.files.classify_summary import ClassifySummaryFileBAC120, ClassifySummaryFileRow

N = 'GCF_000000012.1'
N_TAXONOMY = ['d__Bacteria', 'p__P', 'c__C', 'o__O', 'f__F', 'g__Genus', 's__Genus species']


def make_classify():
    classify = Classify.__new__(Classify)
    classify.logger = logging.getLogger('timestamp')
    classify.gtdb_taxonomy = {'RS_' + N: N_TAXONOMY}
    classify.de_novo_reps_written = set()
    return classify


class TestDeNovoSummaryRows(unittest.TestCase):

    def setUp(self):
        self.dir = tempfile.mkdtemp()
        self.summary = ClassifySummaryFileBAC120(self.dir, 'test')

    def tearDown(self):
        shutil.rmtree(self.dir)

    def add(self, sp_cluster, changing, marker_set_id='bac120', classify=None):
        classify = classify or make_classify()
        return classify._add_de_novo_species_to_summary(
            self.summary, sp_cluster, changing, {}, {}, {}, {}, [], marker_set_id)

    def test_placed_rep_gets_species_and_members_inherit(self):
        row = ClassifySummaryFileRow()
        row.gid = 'U4'
        row.classification = 'd__Bacteria;p__P;c__C;o__O;f__F;g__G2;s__'
        row.note = 'classification based on placement in class-level tree'
        self.summary.add_row(row)
        sp_cluster = {'U4': {'increment': 1, 'radius': 95.5,
                             'members': [('U5', {'ani': 96.2, 'af': 0.9})]}}
        self.add(sp_cluster, {'U5': (N, 'U4')})

        rep = self.summary.get_row('U4')
        self.assertEqual(rep.classification, 'd__Bacteria;p__P;c__C;o__O;f__F;g__G2;s__NovelCluster sp1')
        self.assertIn('Selected to be representative genome for the species cluster', rep.note)
        member = self.summary.get_row('U5')
        self.assertEqual(member.classification, rep.classification)
        self.assertEqual(member.classification_method, 'dereplication_ani')
        self.assertIn('Represented by genome U4 (ANI:96.2%,AF:0.9)', member.note)
        self.assertIn(f'Previously assigned to GTDB representative {N}', member.note)

    def test_rep_moved_from_gtdb_gets_genus_of_previous_species(self):
        sp_cluster = {'Y': {'increment': 2, 'radius': 96.5, 'members': []}}
        self.add(sp_cluster, {'Y': (N, 'Y')})
        row = self.summary.get_row('Y')
        self.assertEqual(row.classification, 'd__Bacteria;p__P;c__C;o__O;f__F;g__Genus;s__NovelCluster sp2')
        self.assertIn(f'Genus from ANI to GTDB representative {N}', row.note)

    def test_genome_moved_from_gtdb_without_species(self):
        self.add({}, {'Y': (N, None)})
        row = self.summary.get_row('Y')
        self.assertEqual(row.classification, 'd__Bacteria;p__P;c__C;o__O;f__F;g__Genus;s__')

    def placed_row(self, gid, classification='d__Bacteria;p__P;c__C;o__O;f__F;g__G2;s__'):
        row = ClassifySummaryFileRow()
        row.gid = gid
        row.classification = classification
        self.summary.add_row(row)

    def test_placed_rep_moved_from_gtdb_is_not_duplicated(self):
        # standalone classify / classify_wf: Y is aligned and placed, so it already has a row
        self.placed_row('Y')
        sp_cluster = {'Y': {'increment': 2, 'radius': 96.5, 'members': []}}
        self.add(sp_cluster, {'Y': (N, 'Y')})
        self.assertEqual(self.summary.get_row('Y').classification,
                         'd__Bacteria;p__P;c__C;o__O;f__F;g__G2;s__NovelCluster sp2')

    def test_placed_genome_without_species_is_not_duplicated(self):
        self.placed_row('Y')
        self.add({}, {'Y': (N, None)})
        self.assertEqual(self.summary.get_row('Y').classification, 'd__Bacteria;p__P;c__C;o__O;f__F;g__G2;s__')

    def test_cluster_of_other_domain_is_skipped(self):
        classify = make_classify()
        sp_cluster = {'R': {'increment': 1, 'radius': 95.0, 'members': [('M', {'ani': 97.0, 'af': 0.8})]}}
        self.add(sp_cluster, {}, marker_set_id='ar53', classify=classify)
        self.assertFalse(self.summary.has_row())
        self.assertEqual(classify.de_novo_reps_written, set())

    def test_unplaced_rep_moved_from_gtdb_only_in_its_domain(self):
        sp_cluster = {'Y': {'increment': 1, 'radius': 96.5, 'members': []}}
        self.add(sp_cluster, {'Y': (N, 'Y')}, marker_set_id='ar53')
        self.assertFalse(self.summary.has_row())

class TestSkaniTriangle(unittest.TestCase):

    def test_run_triangle_parses_sparse_output(self):
        skani = SkANI.__new__(SkANI)
        skani.cpus = 2
        skani.logger = logging.getLogger('timestamp')
        genomes = {'A': '/g/A.fna', 'B': '/g/B.fna', 'C': '/g/C.fna.gz'}
        stdout = ('Ref_file\tQuery_file\tANI\tAlign_fraction_ref\tAlign_fraction_query\tRef_name\tQuery_name\n'
                  '/g/A.fna\t/g/B.fna\t96.80\t99.68\t89.45\tA_c1\tB_c1\n'
                  '/g/A.fna\t/g/C.fna.gz\t95.10\t40.00\t60.00\tA_c1\tC_c1\n')
        proc = mock.Mock(returncode=0, stdout=stdout, stderr='')
        with mock.patch('gtdbtk.external.skani.subprocess.run', return_value=proc) as run:
            results = skani.run_triangle(genomes, 'test')
        args = run.call_args[0][0]
        self.assertEqual(args[:3], ['skani', 'triangle', '-E'])
        self.assertEqual(results, {'B': {'A': {'ani': 96.8, 'af': 0.9968}},
                                   'C': {'A': {'ani': 95.1, 'af': 0.6}}})


if __name__ == '__main__':
    unittest.main()
