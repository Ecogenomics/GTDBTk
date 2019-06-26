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


import unittest


from gtdbtk.biolib_lite.prodigal_biolib import ProdigalGeneFeatureParser


class TestProdigalGeneFeatureParser(unittest.TestCase):

    def setUp(self):
        gff_str = '\n'.join(['##gff-version  3',
                             '# Sequence Data: seqnum=1;seqlen=2001;seqhdr="contig_11394"',
                             '# Model Data: version=Prodigal.v2.6.3;run_type=Single;model="Ab initio";gc_cont=54.92;transl_table=11;uses_sd=1',
                             'contig_11394	Prodigal_v2.6.3	CDS	2	1021	193.2	+	0	ID=1_1;partial=10;start_type=Edge;rbs_motif=None;rbs_spacer=None;gc_cont=0.597;conf=99.99;score=193.23;cscore=190.01;sscore=3.22;rscore=0.00;uscore=0.00;tscore=3.22;',
                             'contig_11394	Prodigal_v2.6.3	CDS	1118	1828	112.4	+	0	ID=1_2;partial=00;start_type=ATG;rbs_motif=None;rbs_spacer=None;gc_cont=0.598;conf=100.00;score=112.42;cscore=110.03;sscore=2.39;rscore=-2.04;uscore=1.59;tscore=2.84;',
                             '# Sequence Data: seqnum=2;seqlen=3000;seqhdr="contig_17797"',
                             '# Model Data: version=Prodigal.v2.6.3;run_type=Single;model="Ab initio";gc_cont=54.92;transl_table=11;uses_sd=1',
                             'contig_17797	Prodigal_v2.6.3	CDS	3	86	3.4	-	0	ID=2_1;partial=10;start_type=ATG;rbs_motif=None;rbs_spacer=None;gc_cont=0.464;conf=68.45;score=3.37;cscore=1.16;sscore=2.21;rscore=-2.04;uscore=2.06;tscore=2.84;',
                             'contig_17797	Prodigal_v2.6.3	CDS	424	1251	136.1	-	0	ID=2_2;partial=00;start_type=ATG;rbs_motif=None;rbs_spacer=None;gc_cont=0.565;conf=100.00;score=136.08;cscore=136.05;sscore=0.03;rscore=-2.04;uscore=-0.12;tscore=2.84;',
                             'contig_17797	Prodigal_v2.6.3	CDS	1384	2811	240.5	+	0	ID=2_3;partial=00;start_type=ATG;rbs_motif=None;rbs_spacer=None;gc_cont=0.551;conf=99.99;score=240.52;cscore=240.99;sscore=-0.47;rscore=-2.04;uscore=-1.27;tscore=2.84;',
                             '# Sequence Data: seqnum=3;seqlen=10;seqhdr="contig_5089"',
                             '# Model Data: version=Prodigal.v2.6.3;run_type=Single;model="Ab initio";gc_cont=44.21;transl_table=4;uses_sd=0',
                             'contig_5089	Prodigal_v2.6.3	CDS	1	4	158.7	+	0	ID=3_1;partial=10;start_type=Edge;rbs_motif=None;rbs_spacer=None;gc_cont=0.503;conf=100.00;score=158.74;cscore=155.52;sscore=3.22;rscore=0.00;uscore=0.00;tscore=3.22;',
                             'contig_5089	Prodigal_v2.6.3	CDS	2	5	161.1	-	0	ID=3_2;partial=00;start_type=ATG;rbs_motif=GGA/GAG/AGG;rbs_spacer=5-10bp;gc_cont=0.527;conf=100.00;score=161.08;cscore=156.63;sscore=4.45;rscore=-0.02;uscore=2.28;tscore=2.84;',
                             'contig_5089	Prodigal_v2.6.3	CDS	7	10	229.5	-	0	ID=3_3;partial=00;start_type=ATG;rbs_motif=None;rbs_spacer=None;gc_cont=0.515;conf=99.99;score=229.54;cscore=231.96;sscore=-2.42;rscore=-2.04;uscore=-4.45;tscore=2.84;',
                             ''])
        self.pgff = ProdigalGeneFeatureParser(gff_str)

    def tearDown(self):
        del self.pgff

    def test__n_seaquences_processed(self):
        """Test that the number of sequences read are correct"""
        self.assertEqual(3, self.pgff.n_sequences_processed())

    def test__build_coding_mask(self):
        """Test that the number of sequences read are correct"""
        expected_mask = [1, 1, 1, 1, 1, 0, 1, 1, 1, 1]
        self.assertEqual(expected_mask, self.pgff.coding_masks['contig_5089'].tolist())

    def test_coding_bases(self):
        """Test that the number of coding bases are calculated correctly"""
        self.assertEqual(0, self.pgff.coding_bases('missing'))
        self.assertEqual(9, self.pgff.coding_bases('contig_5089'))
        self.assertEqual(4, self.pgff.coding_bases('contig_11394', start=2, end=5))
        self.assertEqual(1731, self.pgff.coding_bases('contig_11394'))

    def test_parse_sequence_data(self):
        """Test that the gene data are parsed correctly"""
        self.assertTupleEqual(tuple(self.pgff.sequence_data['contig_11394']), (1, 2001, 'contig_11394'))
        self.assertTupleEqual(tuple(self.pgff.sequence_data['contig_17797']), (2, 3000, 'contig_17797'))
        self.assertTupleEqual(tuple(self.pgff.sequence_data['contig_5089']), (3, 10, 'contig_5089'))

    def test_parse_model_data(self):
        """Test that the gene data are parsed correctly"""
        self.assertTupleEqual(tuple(self.pgff.model_data['contig_11394']),
                              ('Prodigal.v2.6.3', 'Single', 'Ab initio', 54.92, 11, 1))
        self.assertTupleEqual(tuple(self.pgff.model_data['contig_17797']),
                              ('Prodigal.v2.6.3', 'Single', 'Ab initio', 54.92, 11, 1))
        self.assertTupleEqual(tuple(self.pgff.model_data['contig_5089']),
                              ('Prodigal.v2.6.3', 'Single', 'Ab initio', 44.21, 4, 0))

    def test_parse_gene_data(self):
        """Test that the gene data are parsed correctly"""
        self.assertTupleEqual(tuple(self.pgff.called_genes['contig_11394'][0]),
                              ('contig_11394', 'Prodigal_v2.6.3', 'CDS', 2, 1021, 193.23, '+', '0', '1_1', '10', 'Edge',
                               'None', 'None', 0.597, 99.99, 190.01, 3.22, 0.0, 0.0, 3.22))
        self.assertTupleEqual(tuple(self.pgff.called_genes['contig_17797'][1]),
                              (
                                  'contig_17797', 'Prodigal_v2.6.3', 'CDS', 424, 1251, 136.08, '-', '0', '2_2', '00',
                                  'ATG',
                                  'None', 'None', 0.565, 100.0, 136.05, 0.03, -2.04, -0.12, 2.84))
        self.assertTupleEqual(tuple(self.pgff.called_genes['contig_5089'][1]),
                              ('contig_5089', 'Prodigal_v2.6.3', 'CDS', 2, 5, 161.08, '-', '0', '3_2', '00', 'ATG',
                               'GGA/GAG/AGG', '5-10bp', 0.527, 100.0, 156.63, 4.45, -0.02, 2.28, 2.84))

    def test_generate_statistics(self):
        """Test that the statistics are generated correctly"""
        stats = self.pgff.generate_statistics()
        self.assertEqual(round(stats.coding_density, 2), round(4076.0 / (2011 + 3000 + 10), 2))
        self.assertEqual(stats.conf_50, 99.995)
        self.assertEqual(stats.conf_max, 100.0)
        self.assertEqual(round(stats.conf_std, 1), 10.4)
        self.assertEqual(stats.cscore_50, 156.075)
        self.assertEqual(stats.cscore_max, 240.99)
        self.assertEqual(round(stats.cscore_std, 1), 71.2)
        self.assertEqual(stats.gc_cont_50, 0.539)
        self.assertEqual(round(stats.gc_cont_std, 1), 0)
        self.assertEqual(stats.gc_cont_max, 0.598)
        self.assertEqual(stats.genes_called, 8)
        self.assertEqual(stats.rscore_50, -2.04)
        self.assertEqual(round(stats.rscore_std, 1), 1)
        self.assertEqual(stats.rscore_max, 0)
        self.assertEqual(stats.tscore_50, 2.84)
        self.assertEqual(round(stats.tscore_std, 1), 0.2)
        self.assertEqual(stats.tscore_max, 3.22)
        self.assertEqual(stats.uscore_50, 0)
        self.assertEqual(round(stats.uscore_std, 1), 2)
        self.assertEqual(stats.uscore_max, 2.28)


if __name__ == '__main__':
    unittest.main()
