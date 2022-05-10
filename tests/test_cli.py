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

import argparse
import os
import random
import shutil
import string
import tempfile
import unittest

from gtdbtk.biolib_lite.logger import logger_setup
from gtdbtk.config.output import *
from gtdbtk.main import OptionsParser


class TestCli(unittest.TestCase):

    def setUp(self):
        self.identify_dir_reference = os.path.join(os.path.dirname(__file__), 'data/identify_dir_reference/')
        self.align_dir_reference = 'tests/data/align_dir_reference/'
        self.genome_dir = 'gtdbtk/tests/data/genomes/'

        self.options = argparse.ArgumentParser()
        self.options.batchfile = None
        self.options.prefix = 'gtdbtk'
        self.options.cpus = 1
        self.options.extension = 'fna'
        self.options.debug = False
        self.options.force = False
        self.options.genes = False
        self.options.write_single_copy_genes = False

        # align option
        self.options.skip_gtdb_refs = False
        self.options.taxa_filter = None
        self.options.custom_msa_filters = False
        self.options.skip_trimming = False
        self.options.min_consensus = None
        self.options.min_perc_taxa = None
        self.options.skip_gtdb_refs = False
        self.options.cols_per_gene = None
        self.options.max_consensus = None
        self.options.min_perc_aa = 50
        self.options.rnd_seed = 42
        self.options.outgroup_taxon = None

        # classify options
        self.options.scratch_dir = None
        self.options.keep_ref_red = None
        self.options.pplacer_cpus = None
        self.options.min_af = None

        # infer options
        self.options.prot_model = 'WAG'
        self.options.no_support = False
        self.options.no_gamma = True

        self.version = ' unittest'
        self.optionparser = OptionsParser(self.version)
        logger_setup(None, "gtdbtk.log", "GTDB-Tk", self.version, True)
        # self.generic_out_path = 'tests/data/results'
        self.generic_out_path = tempfile.mkdtemp(prefix='gtdbtk_tmp_')

    def tearDown(self):
        shutil.rmtree(self.generic_out_path)


    def test_identify(self):
        tmp_folder = ''.join(random.choice(
            string.ascii_uppercase + string.digits) for _ in range(10))
        identify_options = self.options
        identify_options.genome_dir = self.genome_dir
        identify_options.out_dir = os.path.join(
            self.generic_out_path, tmp_folder, 'identify')
        self.optionparser.identify(identify_options)

        ar53_marker_path = os.path.join(self.options.out_dir,
                                         PATH_AR53_MARKER_SUMMARY.format(prefix=self.options.prefix))

        self.assertTrue(os.path.isfile(
            os.path.join(self.options.out_dir, PATH_BAC120_MARKER_SUMMARY.format(prefix=self.options.prefix))))
        self.assertTrue(os.path.isfile(ar53_marker_path))

        results = {}
        with open(ar53_marker_path, 'r') as f:
            f.readline()
            for line in f:
                infos = line.split('\t', 1)
                results[infos[0]] = infos[1]
        self.assertTrue(results.get('genome_1').startswith('120\t2\t0\t'))

    def test_align(self):
        tmp_folder = ''.join(random.choice(
            string.ascii_uppercase + string.digits) for _ in range(10))
        align_options = self.options
        align_options.identify_dir = self.identify_dir_reference
        align_options.out_dir = os.path.join(
            self.generic_out_path, tmp_folder, 'align')
        self.optionparser.align(align_options)
        path_user_msa = os.path.join(align_options.out_dir, PATH_AR53_USER_MSA.format(prefix=align_options.prefix))
        self.assertTrue(os.path.isfile(path_user_msa))
        with open(path_user_msa, 'r') as f:
            lines = f.read().splitlines()
            last_line = lines[-1]
        self.assertTrue(len(last_line) > 4500)
        self.assertTrue(len(last_line) < 5500)
        self.assertTrue('-' in last_line)
        self.assertFalse(any(char.isdigit() for char in last_line))

    def test_identify_align(self):
        tmp_folder = ''.join(random.choice(
            string.ascii_uppercase + string.digits) for _ in range(10))

        identify_options = self.options
        identify_options.genome_dir = self.genome_dir
        identify_options.out_dir = os.path.join(
            self.generic_out_path, tmp_folder, 'identify')
        self.optionparser.identify(identify_options)

        align_options = self.options
        align_options.identify_dir = os.path.join(
            self.generic_out_path, tmp_folder, 'identify')
        align_options.out_dir = os.path.join(
            self.generic_out_path, tmp_folder, 'align')
        self.optionparser.align(align_options)
        path_user_msa = os.path.join(align_options.out_dir, PATH_AR53_USER_MSA.format(prefix=align_options.prefix))
        self.assertTrue(os.path.isfile(path_user_msa))
        with open(path_user_msa, 'r') as f:
            lines = f.read().splitlines()
            last_line = lines[-1]
        self.assertTrue(len(last_line) > 4500)
        self.assertTrue(len(last_line) < 5500)
        self.assertTrue('-' in last_line)
        self.assertFalse(any(char.isdigit() for char in last_line))

    def test_identify_align_classify(self):
        tmp_folder = ''.join(random.choice(
            string.ascii_uppercase + string.digits) for _ in range(10))

        identify_options = self.options
        identify_options.genome_dir = self.genome_dir
        identify_options.out_dir = os.path.join(
            self.generic_out_path, tmp_folder, 'identify')
        self.optionparser.identify(identify_options)

        align_options = self.options
        align_options.identify_dir = os.path.join(
            self.generic_out_path, tmp_folder, 'identify')
        align_options.out_dir = os.path.join(
            self.generic_out_path, tmp_folder, 'align')
        self.optionparser.align(align_options)
        path_user_msa = os.path.join(align_options.out_dir, PATH_AR53_USER_MSA.format(prefix=align_options.prefix))
        self.assertTrue(os.path.isfile(path_user_msa))
        with open(path_user_msa, 'r') as f:
            lines = f.read().splitlines()
            last_line = lines[-1]
        self.assertTrue(len(last_line) > 4500)
        self.assertTrue(len(last_line) < 5500)
        self.assertTrue('-' in last_line)
        self.assertFalse(any(char.isdigit() for char in last_line))

        classify_options = self.options
        classify_options.genome_dir = self.genome_dir
        classify_options.full_tree = True
        classify_options.align_dir = align_options.out_dir
        classify_options.out_dir = os.path.join(
            self.generic_out_path, tmp_folder, 'classify')
        self.optionparser.classify(classify_options)
        summary_out = os.path.join(classify_options.out_dir,
                                   PATH_AR53_SUMMARY_OUT.format(prefix=classify_options.prefix))
        self.assertTrue(summary_out)
        with open(summary_out, 'r') as f:
            lines = f.read().splitlines()
            last_line = lines[-1]
        infos = last_line.split('\t')
        self.assertEqual(len(infos), 20)
        self.assertTrue(infos[1].startswith('d__Archaea'))

        self.assertTrue(os.path.isdir(os.path.join(classify_options.out_dir, DIR_IDENTIFY_INTERMEDIATE)))
        self.assertTrue(os.path.isdir(os.path.join(classify_options.out_dir, DIR_ALIGN_INTERMEDIATE)))
        self.assertTrue(os.path.isdir(os.path.join(classify_options.out_dir, DIR_CLASSIFY_INTERMEDIATE)))
        self.optionparser.remove_intermediate_files(classify_options.out_dir,'classify_wf')
        self.assertFalse(os.path.exists(os.path.join(classify_options.out_dir, DIR_IDENTIFY_INTERMEDIATE)))
        self.assertFalse(os.path.exists(os.path.join(classify_options.out_dir, DIR_ALIGN_INTERMEDIATE)))
        self.assertFalse(os.path.exists(os.path.join(classify_options.out_dir, DIR_CLASSIFY_INTERMEDIATE)))


    def test_classify_wf(self):
        tmp_folder = ''.join(random.choice(
            string.ascii_uppercase + string.digits) for _ in range(10))
        classify_wf_options = self.options
        classify_wf_options.genome_dir = self.genome_dir
        classify_wf_options.out_dir = os.path.join(
            self.generic_out_path, tmp_folder, 'classify_wf')
        self.optionparser.identify(classify_wf_options)
        classify_wf_options.identify_dir = classify_wf_options.out_dir
        classify_wf_options.align_dir = classify_wf_options.out_dir
        classify_wf_options.taxa_filter = None
        classify_wf_options.custom_msa_filters = False
        classify_wf_options.min_consensus = None
        classify_wf_options.min_perc_taxa = None
        classify_wf_options.skip_gtdb_refs = False
        classify_wf_options.cols_per_gene = None
        classify_wf_options.max_consensus = None
        classify_wf_options.full_tree = True
        self.optionparser.align(classify_wf_options)
        self.optionparser.classify(classify_wf_options)
        summary_out = os.path.join(classify_wf_options.out_dir,
                                   PATH_AR53_SUMMARY_OUT.format(prefix=classify_wf_options.prefix))
        self.assertTrue(os.path.isfile(summary_out))
        with open(summary_out, 'r') as f:
            lines = f.read().splitlines()
            last_line = lines[-1]
        infos = last_line.split('\t')
        self.assertEqual(len(infos), 20)
        self.assertTrue(infos[1].startswith('d__Archaea'))

    def test_infer(self):
        tmp_folder = ''.join(random.choice(
            string.ascii_uppercase + string.digits) for _ in range(10))
        infer_options = self.options
        path_user_msa = PATH_AR53_USER_MSA.format(prefix=self.options.prefix)
        infer_options.msa_file = os.path.join(self.align_dir_reference, path_user_msa)
        infer_options.out_dir = os.path.join(self.generic_out_path, tmp_folder, 'infer')
        infer_options.gamma = False
        # if not os.path.isdir(infer_options.out_dir):
        #     os.makedirs(infer_options.out_dir)
        self.optionparser.infer(infer_options)
        with open(os.path.join(infer_options.out_dir, PATH_TREE_LOG.format(prefix=self.options.prefix)), 'r') as f:
            lines = f.read().splitlines()
            last_line = lines[-1]
        self.assertEqual(last_line.strip(), 'TreeCompleted')
        with open(os.path.join(infer_options.out_dir, PATH_UNROOTED_TREE.format(prefix=self.options.prefix)), 'r') as f:
            lines = f.read().splitlines()
            last_line = lines[-1]
        self.assertTrue('genome_1' in last_line)
        self.assertTrue('genome_2' in last_line)
        self.assertTrue('genome_3' in last_line)

    def test_de_novo_wf(self):
        tmp_folder = ''.join(random.choice(
            string.ascii_uppercase + string.digits) for _ in range(10))
        de_novo_wf_options = self.options
        de_novo_wf_options.genome_dir = self.genome_dir
        de_novo_wf_options.suffix = ".ar53"
        de_novo_wf_options.gamma = False
        de_novo_wf_options.out_dir = os.path.join(
            self.generic_out_path, tmp_folder, 'de_novo_wf')
        de_novo_wf_options.identify_dir = de_novo_wf_options.out_dir
        de_novo_wf_options.msa_file = os.path.join(
            de_novo_wf_options.out_dir, de_novo_wf_options.prefix + de_novo_wf_options.suffix + ".user_msa.fasta")
        self.optionparser.identify(de_novo_wf_options)
        self.optionparser.align(de_novo_wf_options)
        self.optionparser.infer(de_novo_wf_options)

    def test_root(self):
        """Test that rooting is successful when called through the CLI"""
        options = argparse.ArgumentParser()
        options.input_tree = 'tests/data/pplacer_dir_reference/gtdbtk.ar53.classify.tree'
        options.outgroup_taxon = 'p__Altiarchaeota'
        options.output_tree = os.path.join(self.generic_out_path, 'test.rooted.tree')
        options.custom_taxonomy_file = None
        options.gtdbtk_classification_file = None
        self.optionparser.root(options)
        self.assertTrue(os.path.isfile(options.output_tree))


if __name__ == '__main__':
    unittest.main()
