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
import random
import shutil
import string
import unittest

import dendropy

import gtdbtk.config.config as Config
from gtdbtk.biolib_lite.taxonomy import Taxonomy
from gtdbtk.classify import Classify
from gtdbtk.config.output import *


class TestClassify(unittest.TestCase):

    def setUp(self):
        self.classify = Classify()

        self.generic_out_path = 'tests/data/results'
        tmp_folder = ''.join(random.choice(
            string.ascii_uppercase + string.digits) for _ in range(10))
        self.out_dir = os.path.join(self.generic_out_path, tmp_folder)
        if not os.path.exists(self.generic_out_path):
            os.makedirs(self.generic_out_path)
        self.prefix = 'gtdbtk'
        self.pplacer_dir_reference = 'tests/data/pplacer_dir_reference'
        self.aln_dir_ref = 'tests/data/align_dir_reference/align'
        self.user_msa_file = os.path.join(self.aln_dir_ref, 'gtdbtk.ar122.user_msa.fasta')
        self.taxonomy_file = Config.TAXONOMY_FILE
        self.gtdb_taxonomy = Taxonomy().read(self.taxonomy_file)

    def test_standardise_taxonomy(self):
        taxstring = 'p__phylum1;c_class1'
        marker_set = 'bac120'
        new_taxstring = self.classify.standardise_taxonomy(
            taxstring, marker_set)
        self.assertEqual(
            new_taxstring, 'd__Bacteria;p__phylum1;c_class1;o__;f__;g__;s__')

    def test_write_red_dict(self):
        if not os.path.exists(self.out_dir):
            os.makedirs(self.out_dir)
        marker_dict = self.classify._write_red_dict(
            self.out_dir, self.prefix, 'bac120')
        self.assertTrue(len(marker_dict) == 6)
        self.assertTrue('d__' in marker_dict)
        self.assertTrue(marker_dict.get('d__') == 0)
        self.assertTrue('p__' in marker_dict)
        self.assertTrue('c__' in marker_dict)
        self.assertTrue('o__' in marker_dict)
        self.assertTrue('f__' in marker_dict)
        self.assertTrue('g__' in marker_dict)

    def test_get_pplacer_taxonomy(self):
        if not os.path.exists(self.out_dir):
            os.makedirs(self.out_dir)
        tree = dendropy.Tree.get_from_path(os.path.join(os.getcwd(), self.pplacer_dir_reference,
                                                        'gtdbtk.ar122.classify.tree'),
                                           schema='newick',
                                           rooting='force-rooted',
                                           preserve_underscores=True)
        self.classify._get_pplacer_taxonomy(
            self.out_dir, self.prefix, 'ar122', self.user_msa_file, tree)
        results = {}

        with open(os.path.join(self.out_dir, PATH_AR122_PPLACER_CLASS.format(prefix=self.prefix)), 'r') as f:
            for line in f:
                infos = line.strip().split('\t')
                results[infos[0]] = infos[1]
        self.assertTrue(len(results) == 3)
        self.assertTrue('genome_1' in results)
        self.assertTrue('genome_2' in results)
        self.assertTrue('genome_3' in results)
        self.assertEqual(results.get(
            'genome_1'), 'd__Archaea;p__Euryarchaeota;c__Methanobacteria;o__Methanobacteriales;f__Methanobacteriaceae;g__Methanobrevibacter;s__')

    def test_place_genomes(self):
        if not os.path.exists(self.out_dir):
            os.makedirs(self.out_dir)
        tree_file = self.classify.place_genomes(
            self.user_msa_file, 'ar122', self.out_dir, self.prefix)
        with open(tree_file, 'r') as f:
            lines = f.read().splitlines()
            last_line = lines[-1]
        self.assertTrue(last_line.startswith('('))
        self.assertTrue(last_line.endswith('d__Archaea;'))

    def test_formatnote(self):
        first3genomes = list(self.gtdb_taxonomy.keys())[:3]
        sorted_dict = ((first3genomes[0], {'ani': 98.5, 'af': 1.0}), (first3genomes[1], {
            'ani': 92.6, 'af': 1.0}), (first3genomes[2], {'ani': 90.3, 'af': 1.3}))
        labels = [first3genomes[0]]
        note_list = self.classify._formatnote(sorted_dict, labels)
        self.assertTrue(first3genomes[1] in note_list[0])
        self.assertTrue(first3genomes[2] in note_list[1])
        self.assertTrue(note_list[0].endswith(', 92.6, 1.0'))
        self.assertTrue(note_list[1].endswith(', 90.3, 1.3'))

    def test_calculate_red_distances(self):
        tree = os.path.join(self.pplacer_dir_reference,
                            'gtdbtk.ar122.classify.tree')
        result_tree = self.classify._calculate_red_distances(
            tree, self.out_dir)
        egs2 = [eg.length for eg in result_tree.postorder_edge_iter()
                if eg.length is not None]

        self.assertTrue(sum(egs2) / len(egs2) < 0.1)

    def tearDown(self):
        shutil.rmtree(self.generic_out_path)


if __name__ == '__main__':
    unittest.main()
