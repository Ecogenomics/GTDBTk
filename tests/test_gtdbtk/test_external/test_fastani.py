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
import json
from gtdbtk.external.fastani import FastANI
import gtdbtk.config.config as Config

class TestFastANI(unittest.TestCase):

    def setUp(self):
        self.dir_tmp = tempfile.mkdtemp(prefix='gtdbtk_tmp_')
        self.cpus = 1
        self.genome_root = Config.FASTANI_GENOMES

    def tearDown(self):
        shutil.rmtree(self.dir_tmp)

    def test_run(self):
        """Test that FastANI produces the expected output (version dependent)"""

        fa = FastANI(self.cpus, force_single=True)

        """
        a = d__Archaea; p__Altiarchaeota; c__Altiarchaeia; o__Altiarchaeales; f__Altiarchaeaceae; g__Altiarchaeum; s__Altiarchaeum sp001873845
        b = d__Archaea; p__Altiarchaeota; c__Altiarchaeia; o__Altiarchaeales; f__Altiarchaeaceae; g__Altiarchaeum; s__Altiarchaeum sp002083985
        c = d__Bacteria; p__Aquificota; c__Desulfurobacteriia; o__Desulfurobacteriales; f__Desulfurobacteriaceae; g__Thermovibrio; s__Thermovibrio ammonificans
        x = d__Archaea; p__Altiarchaeota; c__Altiarchaeia; o__GCA-002841105; f__GCA-002841105; g__GCA-002841105; s__GCA-002841105 sp002841105
        y = d__Bacteria; p__Aerophobota; c__Aerophobia; o__Aerophobales; f__Aerophobaceae; g__Aerophobus; s__Aerophobus profundus
        z = d__Archaea; p__Halobacterota; c__Archaeoglobi; o__JdFR-21; f__JdFR-21; g__JdFR-21; s__JdFR-21 sp002011165
        """
        d_compare = {'a': {'x', 'y'},
                     'b': {'x'},
                     'c': {'z'}}
        d_paths = {'a': os.path.join(self.genome_root,'GCA/001/873/845', 'GCA_001873845.1_genomic.fna.gz'),
                   'b': os.path.join(self.genome_root,'GCA/002/083/985', 'GCA_002083985.1_genomic.fna.gz'),
                   'c': os.path.join(self.genome_root,'GCF/000/185/805', 'GCF_000185805.1_genomic.fna.gz'),
                   'x': os.path.join(self.genome_root,'GCA/002/841/105', 'GCA_002841105.1_genomic.fna.gz'),
                   'y': os.path.join(self.genome_root,'GCA/000/402/295', 'GCA_000402295.1_genomic.fna.gz'),
                   'z': os.path.join(self.genome_root,'GCA/002/011/165', 'GCA_002011165.1_genomic.fna.gz')}

        result = fa.run(d_compare, d_paths)

        expected = {'a': {'x': {'ani': 82.5201, 'af': 0.57},
                          'y': {'ani': 74.5154, 'af': 0.0}},
                    'b': {'x': {'ani': 84.5846, 'af': 0.55}},
                    'c': {'z': {'ani': 74.4978, 'af': 0.01}}}
        self.assertEqual(json.dumps(result, sort_keys=True), json.dumps(expected, sort_keys=True))

    def test_parse_output_file(self):
        fa = FastANI(self.cpus, force_single=True)

        out_txt = [['q1', 'r1', '83.1234', '5', '10'],
                   ['q1', 'r2', '99.1111', '3', '10'],
                   ['q2', 'r1', '1.12', '555', '1111']]

        expected = {'q1': {'r1': (83.1234, 0.5), 'r2': (99.1111, 0.3)}, 'q2': {'r1': (1.12, 0.5)}}

        path_f1 = os.path.join(self.dir_tmp, 'f1.txt')
        path_f2 = os.path.join(self.dir_tmp, 'f2.txt')
        path_f3 = os.path.join(self.dir_tmp, 'f3.txt')

        with open(path_f1, 'w') as fh:
            for x in out_txt:
                fh.write(' '.join(x) + '\n')

        result_f1 = fa.parse_output_file(path_f1)
        self.assertEqual(result_f1, expected)

        with open(path_f2, 'w') as fh:
            for x in out_txt:
                fh.write('\t'.join(x) + '\n')

        result_f2 = fa.parse_output_file(path_f2)
        self.assertEqual(result_f2, expected)

        open(path_f3, 'w').close()
        result_f3 = fa.parse_output_file(path_f3)
        self.assertEqual(result_f3, {})
