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

import dendropy

import gtdbtk.config.config as Config
from gtdbtk.biolib_lite.taxonomy import Taxonomy
from gtdbtk.classify import Classify
from gtdbtk.config.output import *
from gtdbtk.files.pplacer_classification import PplacerClassifyFileAR53


class TestClassify(unittest.TestCase):

    def setUp(self):
        self.classify = Classify()
        self.out_dir = tempfile.mkdtemp(prefix='gtdbtk_tmp_')
        self.prefix = 'gtdbtk'
        self.pplacer_dir_reference = 'tests/data/pplacer_dir_reference'
        self.aln_dir_ref = 'tests/data/align_dir_reference/align'
        self.user_msa_file = os.path.join(self.aln_dir_ref, 'gtdbtk.ar53.user_msa.fasta')
        self.taxonomy_file = Config.TAXONOMY_FILE
        self.gtdb_taxonomy = Taxonomy().read(self.taxonomy_file)

    def tearDown(self):
        shutil.rmtree(self.out_dir)

    def test_standardise_taxonomy(self):
        taxstring = 'p__phylum1;c_class1'
        marker_set = 'bac120'
        new_taxstring = self.classify.standardise_taxonomy(
            taxstring, marker_set)
        self.assertEqual(
            new_taxstring, 'd__Bacteria;p__phylum1;c_class1;o__;f__;g__;s__')

        # Test that the correct domain is returned.
        self.assertEqual(self.classify.standardise_taxonomy('p__P;c__C;o__O;f__F;g__G;s__S', 'bac120'),
                         'd__Bacteria;p__P;c__C;o__O;f__F;g__G;s__S')
        self.assertEqual(self.classify.standardise_taxonomy('p__P;c__C;o__O;f__F;g__G;s__S', 'ar53'),
                         'd__Archaea;p__P;c__C;o__O;f__F;g__G;s__S')

        # Remove ranks and check
        rank_order = {'p': 0, 'c': 1, 'o': 2, 'f': 3, 'g': 4, 's': 5}
        rank_lst = ['p__P', 'c__C', 'o__O', 'f__F', 'g__G', 's__S']
        ranks = {'p': 'P', 'c': 'C', 'o': 'O', 'f': 'F', 'g': 'G', 's': 'S'}
        dom_info = {'d__Bacteria': 'bac120', 'd__Archaea': 'ar53'}

        for k in range(1, len(ranks) - 1):
            for cur_domain in ('d__Bacteria', 'd__Archaea'):
                ranks_selected = rank_lst[0:-k]
                expected = list()
                test_lst = list()
                for cur_rank, _ in sorted(rank_order.items(), key=lambda x: [1]):
                    if cur_rank in ranks_selected:
                        test_lst.append(f'{cur_rank}__{ranks[cur_rank]}')
                        expected.append(f'{cur_rank}__{ranks[cur_rank]}')
                    else:
                        expected.append(f'{cur_rank}__')

                expected_str = f'{cur_domain};{";".join(expected)}'
                test_str = ";".join(test_lst)

                cur_dom = dom_info[cur_domain]
                test_value = self.classify.standardise_taxonomy(test_str, cur_dom)
                self.assertEqual(expected_str, test_value)

    def test_get_pplacer_taxonomy(self):
        if not os.path.exists(self.out_dir):
            os.makedirs(self.out_dir)
        tree = dendropy.Tree.get_from_path(os.path.join(os.getcwd(), self.pplacer_dir_reference,
                                                        'gtdbtk.ar53.classify.tree'),
                                           schema='newick',
                                           rooting='force-rooted',
                                           preserve_underscores=True)
        pplacer_classify_file = PplacerClassifyFileAR53(self.out_dir, self.prefix)
        self.classify._get_pplacer_taxonomy(pplacer_classify_file, 'ar53', self.user_msa_file, tree)
        results = {}

        with open(os.path.join(self.out_dir, PATH_AR53_PPLACER_CLASS.format(prefix=self.prefix)), 'r') as f:
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
        tree_file = self.classify.place_genomes(
            self.user_msa_file, 'ar53', self.out_dir, self.prefix)
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
                            'gtdbtk.ar53.classify.tree')
        result_tree = self.classify._calculate_red_distances(
            tree, self.out_dir)
        egs2 = [eg.length for eg in result_tree.postorder_edge_iter()
                if eg.length is not None]

        self.assertTrue(sum(egs2) / len(egs2) < 0.1)


if __name__ == '__main__':
    unittest.main()
