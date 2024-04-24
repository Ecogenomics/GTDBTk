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
from gtdbtk.config.common import CONFIG
from gtdbtk.external.skani import SkANI


class TestSkANI(unittest.TestCase):

    def setUp(self):
        self.dir_tmp = tempfile.mkdtemp(prefix='gtdbtk_tmp_')
        self.cpus = 1
        self.genome_root = CONFIG.SKANI_GENOMES

    def tearDown(self):
        shutil.rmtree(self.dir_tmp)

    def test_run(self):
        """Test that skani produces the expected output (version dependent)"""

        fa = SkANI(self.cpus, force_single=True)

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
