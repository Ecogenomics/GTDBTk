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

import gtdbtk.config.config as Config
from gtdbtk.external.fasttree import FastTree


class TestFastTree(unittest.TestCase):

    def setUp(self):
        self.dir_tmp = tempfile.mkdtemp(prefix='gtdbtk_tmp_')
        self.cpus = 1
        self.genome_root = Config.FASTANI_GENOMES
        self.test_msa = 'tests/data/example.faa'
        self.test_msa_gz = 'tests/data/example.faa.gz'

    def tearDown(self):
        shutil.rmtree(self.dir_tmp)

    def test_run(self):
        ft = FastTree()
        output_tree = os.path.join(self.dir_tmp, 'output.tree')
        tree_log = os.path.join(self.dir_tmp, 'tree.log')
        fasttree_log = os.path.join(self.dir_tmp, 'fasttree.log')
        prot_model = 'WAG'
        no_support = False
        gamma = False
        msa_file = self.test_msa
        cpus = self.cpus
        ft.run(output_tree, tree_log, fasttree_log, prot_model, no_support, gamma, msa_file, cpus)
        self.assertTrue(os.path.isfile(output_tree))

    def test_run_gzip(self):
        ft = FastTree()
        output_tree = os.path.join(self.dir_tmp, 'output.tree')
        tree_log = os.path.join(self.dir_tmp, 'tree.log')
        fasttree_log = os.path.join(self.dir_tmp, 'fasttree.log')
        prot_model = 'WAG'
        no_support = False
        gamma = False
        msa_file = self.test_msa_gz
        cpus = self.cpus
        ft.run(output_tree, tree_log, fasttree_log, prot_model, no_support, gamma, msa_file, cpus)
        self.assertTrue(os.path.isfile(output_tree))
