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
import gtdbtk
import argparse
from gtdbtk.main import OptionsParser
from gtdbtk import tools
import shutil
import os
import logging
from gtdbtk.biolib_lite.logger import logger_setup
import string
import random


class TestCli(unittest.TestCase):

    def setUp(self):
        self.identify_dir_reference = 'tests/data/identify_dir_reference/'
        self.align_dir_reference = 'tests/data/align_dir_reference/'
        self.genome_dir = 'gtdbtk/tests/data/genomes/'

        self.options = argparse.ArgumentParser()
        self.options.batchfile = None
        self.options.prefix = 'gtdbtk'
        self.options.cpus = 1
        self.options.extension = 'fna'
        self.options.debug = False

        # align option
        self.options.skip_gtdb_refs = False
        self.options.taxa_filter = None
        self.options.custom_msa_filters = False
        self.options.min_consensus = None
        self.options.min_perc_taxa = None
        self.options.skip_gtdb_refs = False
        self.options.cols_per_gene = None
        self.options.max_consensus = None
        self.options.min_perc_aa = 50

        # classify options
        self.options.scratch_dir = None

        # infer options
        self.options.prot_model = 'WAG'
        self.options.no_support = False
        self.options.no_gamma = True

        self.version = ' unittest'
        self.optionparser = OptionsParser(self.version)
        logger_setup(None, "gtdbtk.log", "GTDB-Tk", self.version, True)
        self.generic_out_path = 'tests/data/results'

    def test_identify(self):
        tmp_folder = ''.join(random.choice(
            string.ascii_uppercase + string.digits) for _ in range(10))
        identify_options = self.options
        identify_options.genome_dir = self.genome_dir
        identify_options.out_dir = os.path.join(
            self.generic_out_path, tmp_folder, 'identify')
        self.optionparser.identify(identify_options)
        self.assertTrue(os.path.isfile(os.path.join(
            self.options.out_dir, 'gtdbtk_bac120_markers_summary.tsv')))
        self.assertTrue(os.path.isfile(os.path.join(
            self.options.out_dir, 'gtdbtk_ar122_markers_summary.tsv')))

        results = {}
        with open(os.path.join(identify_options.out_dir, 'gtdbtk_ar122_markers_summary.tsv'), 'r') as f:
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
        self.assertTrue(os.path.isfile(os.path.join(
            align_options.out_dir, 'gtdbtk.ar122.user_msa.fasta')))
        with open(os.path.join(align_options.out_dir, 'gtdbtk.ar122.user_msa.fasta'), 'r') as f:
            lines = f.read().splitlines()
            last_line = lines[-1]
        self.assertTrue(len(last_line) > 4500)
        self.assertTrue(len(last_line) < 5500)
        self.assertTrue('-' in last_line)
        self.assertFalse(any(char.isdigit() for char in last_line))

    def test_classify(self):
        tmp_folder = ''.join(random.choice(
            string.ascii_uppercase + string.digits) for _ in range(10))
        classify_options = self.options
        classify_options.genome_dir = self.genome_dir
        classify_options.align_dir = self.align_dir_reference
        classify_options.out_dir = os.path.join(
            self.generic_out_path, tmp_folder, 'classify')
        self.optionparser.classify(classify_options)
        self.assertTrue(os.path.isfile(os.path.join(
            classify_options.out_dir, 'gtdbtk.ar122.summary.tsv')))
        with open(os.path.join(classify_options.out_dir, 'gtdbtk.ar122.summary.tsv'), 'r') as f:
            lines = f.read().splitlines()
            last_line = lines[-1]
        infos = last_line.split('\t')
        self.assertEquals(len(infos), 16)
        self.assertTrue(infos[1].startswith('d__Archaea'))

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
        self.assertTrue(os.path.isfile(os.path.join(
            align_options.out_dir, 'gtdbtk.ar122.user_msa.fasta')))
        with open(os.path.join(align_options.out_dir, 'gtdbtk.ar122.user_msa.fasta'), 'r') as f:
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
        self.assertTrue(os.path.isfile(os.path.join(
            align_options.out_dir, 'gtdbtk.ar122.user_msa.fasta')))
        with open(os.path.join(align_options.out_dir, 'gtdbtk.ar122.user_msa.fasta'), 'r') as f:
            lines = f.read().splitlines()
            last_line = lines[-1]
        self.assertTrue(len(last_line) > 4500)
        self.assertTrue(len(last_line) < 5500)
        self.assertTrue('-' in last_line)
        self.assertFalse(any(char.isdigit() for char in last_line))

        classify_options = self.options
        classify_options.genome_dir = self.genome_dir
        classify_options.align_dir = align_options.out_dir
        classify_options.out_dir = os.path.join(
            self.generic_out_path, tmp_folder, 'classify')
        self.optionparser.classify(classify_options)
        self.assertTrue(os.path.isfile(os.path.join(
            classify_options.out_dir, 'gtdbtk.ar122.summary.tsv')))
        with open(os.path.join(classify_options.out_dir, 'gtdbtk.ar122.summary.tsv'), 'r') as f:
            lines = f.read().splitlines()
            last_line = lines[-1]
        infos = last_line.split('\t')
        self.assertEquals(len(infos), 16)
        self.assertTrue(infos[1].startswith('d__Archaea'))

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
        self.optionparser.align(classify_wf_options)
        self.optionparser.classify(classify_wf_options)

        self.assertTrue(os.path.isfile(os.path.join(
            classify_wf_options.out_dir, 'gtdbtk.ar122.summary.tsv')))
        with open(os.path.join(classify_wf_options.out_dir, 'gtdbtk.ar122.summary.tsv'), 'r') as f:
            lines = f.read().splitlines()
            last_line = lines[-1]
        infos = last_line.split('\t')
        self.assertEquals(len(infos), 16)
        self.assertTrue(infos[1].startswith('d__Archaea'))

    def test_infer(self):
        tmp_folder = ''.join(random.choice(
            string.ascii_uppercase + string.digits) for _ in range(10))
        infer_options = self.options
        infer_options.msa_file = os.path.join(
            self.align_dir_reference, 'gtdbtk.ar122.user_msa.fasta')
        infer_options.out_dir = os.path.join(
            self.generic_out_path, tmp_folder, 'infer')
        self.optionparser.infer(infer_options)
        with open(os.path.join(infer_options.out_dir, 'gtdbtk.tree.log'), 'r') as f:
            lines = f.read().splitlines()
            last_line = lines[-1]
        self.assertEqual(last_line.strip(), 'TreeCompleted')
        with open(os.path.join(infer_options.out_dir, 'gtdbtk.unrooted.tree'), 'r') as f:
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
        de_novo_wf_options.suffix = ".ar122"
        de_novo_wf_options.out_dir = os.path.join(
            self.generic_out_path, tmp_folder, 'de_novo_wf')
        de_novo_wf_options.identify_dir = de_novo_wf_options.out_dir
        de_novo_wf_options.msa_file = os.path.join(
            de_novo_wf_options.out_dir, de_novo_wf_options.prefix + de_novo_wf_options.suffix + ".user_msa.fasta")
        self.optionparser.identify(de_novo_wf_options)
        self.optionparser.align(de_novo_wf_options)
        self.optionparser.infer(de_novo_wf_options)

    def tearDown(self):
        shutil.rmtree(self.generic_out_path)


if __name__ == '__main__':
    unittest.main()
