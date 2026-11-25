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

"""Tests for de novo species clustering (--de_novo_species), using the cases
from the design slides. No reference data or third-party binaries are needed.
"""

import gzip
import os
import shutil
import tempfile
import unittest
from unittest import mock

from gtdbtk.biolib_lite.common import canonical_gid
from gtdbtk.de_novo_species import GenomeClustering, QCedGenome
from gtdbtk.exceptions import GTDBTkExit, GenomeNotInQCFile
from gtdbtk.files.classify_summary import ClassifySummaryFileRow

AF = 0.9

# GTDB representatives
S1, S2, S3, S4 = 'GCF_000000001.1', 'GCF_000000002.1', 'GCF_000000003.1', 'GCF_000000004.1'
G, G2, N = 'GCF_000000010.1', 'GCF_000000011.1', 'GCF_000000012.1'


class FakeRadii(object):
    """Stands in for GTDBRadiiFile: radius per GTDB representative."""

    def __init__(self, radii):
        self.radii = {canonical_gid(gid): radius for gid, radius in radii.items()}

    def get_rep_ani(self, rep):
        return self.radii[rep]


def hit(ani, af=AF):
    return {'ani': ani, 'af': af}


class ClusteringTestCase(unittest.TestCase):

    def setUp(self):
        self.dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.dir)

    def make_clustering(self, qc_values, radii, qc_text=None):
        """qc_values: dict[gid] = (completeness, contamination)."""
        genome_path = {}
        for gid in qc_values:
            path = os.path.join(self.dir, f'{gid}.fna')
            with open(path, 'w') as fh:
                fh.write(f'>{gid}_1\nACGTACGTACGT\n')
            genome_path[gid] = path
        if qc_text is None:
            qc_text = 'Name\tCompleteness\tContamination\n' + ''.join(
                f'{gid}\t{comp}\t{cont}\n' for gid, (comp, cont) in qc_values.items())
        qc_path = os.path.join(self.dir, 'qc.tsv')
        with open(qc_path, 'w') as fh:
            fh.write(qc_text)
        clustering = GenomeClustering(qc_path, genome_path, FakeRadii(radii), self.dir, 'test')
        clustering.parse()
        return clustering

    @staticmethod
    def build_hits(pairs):
        ani_results = {}
        for gid_a, gid_b, ani, af in pairs:
            GenomeClustering.add_hit(ani_results, gid_a, gid_b, hit(ani, af))
        return ani_results


class TestSlidesCases(ClusteringTestCase):
    """GTDB S1-S4 and user genomes U1-U5 from the slides."""

    RADII = {S1: 95.0, S2: 95.0, S3: 96.0, S4: 96.0}
    # Step 1 (ANI screen): U1 -> S1, U2 -> S2, U5 -> S4
    GTDB_ASSIGNMENT = {'U1': (S1, hit(96.5)), 'U2': (S2, hit(96.7)), 'U5': (S4, hit(96.1))}

    def run_case(self, u3_completeness, u4_completeness):
        qc_values = {g: (90.0, 1.0) for g in ('U1', 'U2', 'U5')}
        qc_values['U3'] = (u3_completeness, 1.0)
        qc_values['U4'] = (u4_completeness, 1.0)
        clustering = self.make_clustering(qc_values, self.RADII)
        ani_results = self.build_hits([
            ('U5', 'U4', 96.2, AF),
            ('U4', S4, 95.5, AF),
            ('U4', 'U3', 95.4, AF),
            ('U3', S2, 91.0, AF),
        ])
        return clustering.cluster(dict(self.GTDB_ASSIGNMENT), ani_results)

    def test_case_a_u4_scores_higher(self):
        sp_cluster, _, _, _ = self.run_case(80, 95)
        self.assertEqual(set(sp_cluster), {'U4', 'U3'})
        self.assertEqual(sp_cluster['U4']['radius'], 95.5)  # ANI to S4
        self.assertEqual(sp_cluster['U3']['radius'], 95.4)  # ANI to U4 (D1)
        self.assertEqual(sp_cluster['U4']['increment'], 1)
        self.assertEqual(sp_cluster['U3']['increment'], 2)

    def test_case_b_u3_scores_higher(self):
        sp_cluster, _, _, _ = self.run_case(95, 80)
        self.assertEqual(set(sp_cluster), {'U3', 'U4'})
        self.assertEqual(sp_cluster['U3']['radius'], 95.0)  # nothing within 95%
        # U4 is inside U3's 95% (95.4) but its closest is S4 (95.5): 2.8 rule -> new rep
        self.assertEqual(sp_cluster['U4']['radius'], 95.5)

    def test_final_pass_moves_u5_to_d1(self):
        for u3, u4 in ((80, 95), (95, 80)):
            sp_cluster, gtdb_members, changing, unassigned = self.run_case(u3, u4)
            self.assertEqual(sp_cluster['U4']['members'], [('U5', hit(96.2))])
            self.assertEqual(changing, {'U5': (S4, 'U4')})
            self.assertEqual({g: r for g, (r, _h) in gtdb_members.items()}, {'U1': S1, 'U2': S2})
            self.assertEqual(unassigned, set())


class TestFuzzySpecies(ClusteringTestCase):
    """X/Y example: G and G' are fuzzy (98%), N has 95%."""

    RADII = {G: 98.0, G2: 98.0, N: 95.0}

    def run_case(self, y_completeness=90):
        clustering = self.make_clustering({'X': (95, 1), 'Y': (y_completeness, 1)}, self.RADII)
        ani_results = self.build_hits([
            ('X', G, 97.0, AF),
            ('X', N, 96.0, AF),
            ('X', G2, 96.0, AF),
            ('Y', 'X', 96.5, AF),
        ])
        # Step 1 put Y in N (95.5%)
        return clustering.cluster({'Y': (N, hit(95.5))}, ani_results)

    def test_x_radius_is_ani_to_closest(self):
        sp_cluster, _, _, _ = self.run_case()
        self.assertEqual(sp_cluster['X']['radius'], 97.0)  # not G's 98

    def test_y_loses_gtdb_species_and_becomes_rep(self):
        sp_cluster, gtdb_members, changing, unassigned = self.run_case()
        self.assertEqual(sp_cluster['Y']['radius'], 96.5)
        self.assertEqual(sp_cluster['Y']['increment'], 2)  # created in the second round
        self.assertEqual(changing, {'Y': (N, 'Y')})
        self.assertEqual(gtdb_members, {})
        self.assertEqual(unassigned, set())

    def test_y_failing_gate_has_no_species(self):
        sp_cluster, _, changing, unassigned = self.run_case(y_completeness=40)
        self.assertNotIn('Y', sp_cluster)
        self.assertEqual(unassigned, {'Y'})
        self.assertEqual(changing, {'Y': (N, None)})


class TestRules(ClusteringTestCase):

    def test_low_quality_never_representative(self):
        clustering = self.make_clustering({'Z': (40, 1)}, {S1: 95.0})
        sp_cluster, _, _, unassigned = clustering.cluster({}, {})
        self.assertEqual(sp_cluster, {})
        self.assertEqual(unassigned, {'Z'})

    def test_low_quality_joins_cluster(self):
        clustering = self.make_clustering({'A': (95, 1), 'Z': (40, 1)}, {})
        ani_results = self.build_hits([('A', 'Z', 97.0, AF)])
        sp_cluster, _, _, _ = clustering.cluster({}, ani_results)
        self.assertEqual(set(sp_cluster), {'A'})
        self.assertEqual(sp_cluster['A']['members'], [('Z', hit(97.0))])

    def test_contamination_gate(self):
        clustering = self.make_clustering({'A': (95, 12), 'B': (95, 5)}, {})
        self.assertFalse(clustering.qc_genomes['A'].can_be_representative())
        self.assertTrue(clustering.qc_genomes['B'].can_be_representative())

    def test_best_score_becomes_rep(self):
        clustering = self.make_clustering({'A': (80, 1), 'B': (95, 1)}, {})
        ani_results = self.build_hits([('A', 'B', 98.0, AF)])
        sp_cluster, _, _, _ = clustering.cluster({}, ani_results)
        self.assertEqual(set(sp_cluster), {'B'})
        self.assertEqual(sp_cluster['B']['members'], [('A', hit(98.0))])

    def test_low_af_hits_ignored(self):
        clustering = self.make_clustering({'U': (95, 1)}, {S1: 95.0})
        ani_results = self.build_hits([('U', S1, 99.0, 0.3)])
        sp_cluster, _, _, _ = clustering.cluster({}, ani_results)
        self.assertEqual(sp_cluster['U']['radius'], 95.0)

    def test_closest_tie_broken_by_af_then_id(self):
        clustering = self.make_clustering({'U': (95, 1)}, {S1: 95.0, S2: 95.0, S3: 95.0})
        ani_results = self.build_hits([('U', S2, 96.0, 0.8), ('U', S1, 96.0, 0.8), ('U', S3, 96.0, 0.7)])
        self.assertEqual(clustering.closest_representative('U', ani_results, {})[0], S1)

    def test_add_hit_keeps_best_direction(self):
        ani_results = {}
        GenomeClustering.add_hit(ani_results, 'A', 'B', hit(96.0, 0.6))
        GenomeClustering.add_hit(ani_results, 'B', 'A', hit(96.2, 0.6))
        GenomeClustering.add_hit(ani_results, 'A', 'B', hit(96.1, 0.6))
        self.assertEqual(ani_results['A']['B'], hit(96.2, 0.6))
        self.assertEqual(ani_results['B']['A'], hit(96.2, 0.6))

    def test_input_order_does_not_matter(self):
        pairs = [('A', 'B', 96.0, AF), ('B', 'C', 96.0, AF), ('A', 'C', 95.2, AF)]
        results = []
        for ordered in (pairs, list(reversed(pairs))):
            clustering = self.make_clustering({'A': (90, 1), 'B': (90, 1), 'C': (90, 1)}, {})
            sp_cluster, _, _, _ = clustering.cluster({}, self.build_hits(ordered))
            results.append({k: (v['radius'], v['members']) for k, v in sp_cluster.items()})
        self.assertEqual(results[0], results[1])

    def test_cluster_without_parse(self):
        clustering = GenomeClustering(os.path.join(self.dir, 'qc.tsv'), {}, FakeRadii({}), self.dir, 'test')
        with self.assertRaises(GTDBTkExit):
            clustering.cluster({}, {})

    def test_gtdb_assignment_from_rows(self):
        row = ClassifySummaryFileRow()
        row.gid, row.closest_genome_ref = 'U1', S1
        row.closest_genome_ani, row.closest_genome_af = 96.5, 0.9
        self.assertEqual(GenomeClustering.gtdb_assignment_from_rows([row]), {'U1': (S1, hit(96.5, 0.9))})

    def test_compute_ani_uses_skani_triangle(self):
        clustering = self.make_clustering({'A': (90, 1), 'B': (90, 1), 'C': (90, 1)}, {})
        calls = {}

        class FakeSkANI(object):
            def __init__(self, cpus, force_single):
                pass

            def run_triangle(self, genomes, prefix, skani_preset=None):
                calls['genomes'] = genomes
                return {'B': {'A': hit(96.0)}, 'C': {'A': hit(95.5), 'B': hit(97.0)}}

        with mock.patch('gtdbtk.de_novo_species.SkANI', FakeSkANI):
            ani_results = clustering.compute_ani(['C', 'A', 'B'])
        self.assertEqual(sorted(calls['genomes']), ['A', 'B', 'C'])
        self.assertEqual(ani_results['A'], {'B': hit(96.0), 'C': hit(95.5)})
        self.assertEqual(ani_results['C']['B'], hit(97.0))

    def test_compute_ani_single_genome(self):
        clustering = self.make_clustering({'A': (90, 1)}, {})
        self.assertEqual(clustering.compute_ani(['A']), {})

    def test_write_and_read_clusters(self):
        clustering = self.make_clustering({'A': (95, 1)}, {})
        sp_cluster = {'A': {'increment': 1, 'radius': 95.0, 'members': [('Z', hit(97.0))]}}
        changing = {'Z': (S1, 'A'), 'Y': (N, None)}
        path = os.path.join(self.dir, 'clusters.tsv')
        clustering.write_clusters(path, sp_cluster, changing)
        with open(path) as fh:
            lines = [line.rstrip('\n').split('\t') for line in fh]
        self.assertEqual(lines[1], ['A', 'A', '1', 's__NovelCluster sp1', '95.0', 'N/A', 'N/A', 'N/A'])
        self.assertEqual(lines[2], ['Z', 'A', '1', 's__NovelCluster sp1', '95.0', '97.0', '0.9', S1])
        self.assertEqual(lines[3], ['Y', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', N])
        read_cluster, read_changing = GenomeClustering.read_clusters(path)
        self.assertEqual(read_cluster, {'A': {'increment': 1, 'radius': 95.0, 'members': [('Z', hit(97.0))]}})
        self.assertEqual(read_changing, changing)

    def test_gtdb_assignment_from_ani_screen(self):
        classified = {'bac120': {'U1': {S1: {'ani': 96.5, 'af': 0.9, 'other_refs': 'x'}}},
                      'ar53': {'U2': {S2: {'ani': 97.0, 'af': 0.8}}}}
        self.assertEqual(GenomeClustering.gtdb_assignment_from_ani_screen(classified),
                         {'U1': (S1, hit(96.5, 0.9)), 'U2': (S2, hit(97.0, 0.8))})
        self.assertEqual(GenomeClustering.gtdb_assignment_from_ani_screen(None), {})

class TestQCFile(ClusteringTestCase):

    def test_checkm1(self):
        text = 'Bin Id\tMarker lineage\tCompleteness\tContamination\ng1\tk__Bacteria\t97.5\t1.2\n'
        clustering = self.make_clustering({'g1': None}, {}, qc_text=text)
        self.assertEqual(clustering.de_novo_file_format, 1)
        self.assertEqual((clustering.qc_genomes['g1'].completeness, clustering.qc_genomes['g1'].contamination),
                         (97.5, 1.2))

    def test_checkm2_both_models(self):
        text = ('Name\tCompleteness_General\tContamination\tCompleteness_Specific\tCompleteness_Model_Used\n'
                'g1\t90.0\t2.0\t95.0\tNeural Network (Specific Model)\n'
                'g2\t80.0\t3.0\t70.0\tGradient Boost (General Model)\n')
        clustering = self.make_clustering({'g1': None, 'g2': None}, {}, qc_text=text)
        self.assertEqual(clustering.qc_genomes['g1'].completeness, 95.0)
        self.assertEqual(clustering.qc_genomes['g2'].completeness, 80.0)

    def test_checkm2_one_model(self):
        text = 'Name\tCompleteness_General\tContamination\ng1\t88\t4\n'
        clustering = self.make_clustering({'g1': None}, {}, qc_text=text)
        self.assertEqual(clustering.qc_genomes['g1'].completeness, 88.0)

    def test_generic_csv_and_extra_rows(self):
        text = 'Name,Completeness,Contamination\ng1,99,0.5\nnot_a_user_genome,50,5\n\n'
        clustering = self.make_clustering({'g1': None}, {}, qc_text=text)
        self.assertEqual(set(clustering.qc_genomes), {'g1'})

    def test_unknown_format(self):
        with self.assertRaises(GTDBTkExit):
            self.make_clustering({'g1': None}, {}, qc_text='id\tcomp\tcont\ng1\t99\t1\n')

    def test_missing_genome(self):
        with self.assertRaises(GenomeNotInQCFile):
            self.make_clustering({'g1': None, 'g2': None}, {},
                                 qc_text='Name\tCompleteness\tContamination\ng1\t99\t1\n')

    def test_seq_stats_and_score(self):
        path = os.path.join(self.dir, 'g.fna.gz')
        with gzip.open(path, 'wt') as fh:
            fh.write('>s1\nACGT' + 'N' * 10 + 'ACGTNN\n>s2\nNNACGT\n')
        qg = QCedGenome('g', path, 90.0, 2.0)
        self.assertEqual((qg.contig_count, qg.ambiguous_bases), (3, 14))
        self.assertAlmostEqual(qg.score, 90 - 10 - 5 * 3 / 100 - 5 * 14 / 10000)


if __name__ == '__main__':
    unittest.main()
