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
import re
import shutil
import tempfile
import hashlib
import unittest

from gtdbtk.exceptions import *
from gtdbtk.main import OptionsParser
from gtdbtk.biolib_lite.exceptions import *

class TestOptionsParser(unittest.TestCase):

    def setUp(self):
        self.options_parser = OptionsParser(-1)

        pass

    def tearDown(self):
        pass

    def test__verify_genome_id__valid(self):
        """ Test that a valid genome id returns True. """
        self.assertTrue(self.options_parser._verify_genome_id('genome_1'))

    def test__verify_genome_id__invalid(self):
        """ Test that invalid genome ids throw an exception. """
        for c in list('()[],;='):
            self.assertRaises(GenomeNameInvalid, self.options_parser._verify_genome_id, 'genome%s1' % c)

    def test__genomes_to_process__genome_dir__valid(self):
        """ Test that the expected results are returned when using genome_dir. """
        try:
            tmp_genome_dir = tempfile.mkdtemp()
            open(os.path.join(tmp_genome_dir, 'genome_1.fna'), 'a').close()
            open(os.path.join(tmp_genome_dir, 'genome_2.fna'), 'a').close()
            open(os.path.join(tmp_genome_dir, 'other_file.txt'), 'a').close()
            results = self.options_parser._genomes_to_process(tmp_genome_dir, '', 'fna')
            expected = {'genome_1': os.path.join(tmp_genome_dir, 'genome_1.fna'),
                        'genome_2': os.path.join(tmp_genome_dir, 'genome_2.fna')}
            self.assertDictEqual(results, expected)
        finally:
            shutil.rmtree(tmp_genome_dir)

    def test__genomes_to_process__batchfile__valid(self):
        """ Test that the expected results are returned when using batchfile """
        try:
            tmp_genome_dir = tempfile.mkdtemp()
            path_batchfile = os.path.join(tmp_genome_dir, 'batchfile.txt')
            path_genome_1 = os.path.join(tmp_genome_dir, 'genome_1.fna')
            path_genome_2 = os.path.join(tmp_genome_dir, 'genome_2.fna')
            open(path_genome_1, 'a').close()
            open(path_genome_2, 'a').close()

            with open(path_batchfile, 'a') as f:
                f.write('%s\tgenome_1\n' % path_genome_1)
                f.write('\n')
                f.write('%s\tgenome_2\n' % path_genome_2)

            results = self.options_parser._genomes_to_process('', path_batchfile, 'fna')
            expected = {'genome_1': path_genome_1, 'genome_2': path_genome_2}
            self.assertDictEqual(results, expected)
        finally:
            shutil.rmtree(tmp_genome_dir)

    def test__genomes_to_process__batchfile__invalid_columns(self):
        """ Test that a batchfile containing columns not equal to 2 throws an exception. """
        try:
            tmp_genome_dir = tempfile.mkdtemp()
            path_batchfile = os.path.join(tmp_genome_dir, 'batchfile.txt')
            path_genome_1 = os.path.join(tmp_genome_dir, 'genome_1.fna')
            path_genome_2 = os.path.join(tmp_genome_dir, 'genome_2.fna')
            open(path_genome_1, 'a').close()
            open(path_genome_2, 'a').close()

            with open(path_batchfile, 'a') as f:
                f.write('%s\tgenome_1\n' % path_genome_1)
                f.write('\n')
                f.write('%s\tgenome_2\tfoo\n' % path_genome_2)

            self.assertRaises(GenomeBatchfileMalformed, self.options_parser._genomes_to_process, '', path_batchfile,
                              'fna')
        finally:
            shutil.rmtree(tmp_genome_dir)

    def test__genomes_to_process__batchfile__blank_genome_path(self):
        """ Test that a batchfile containing a blank genome path throws an exception. """
        try:
            tmp_genome_dir = tempfile.mkdtemp()
            path_batchfile = os.path.join(tmp_genome_dir, 'batchfile.txt')
            path_genome_1 = os.path.join(tmp_genome_dir, 'genome_1.fna')
            path_genome_2 = os.path.join(tmp_genome_dir, 'genome_2.fna')
            open(path_genome_1, 'a').close()
            open(path_genome_2, 'a').close()

            with open(path_batchfile, 'a') as f:
                f.write('%s\tgenome_1\n' % path_genome_1)
                f.write('\n')
                f.write('%s\tgenome_2\n' % '')

            self.assertRaises(GenomeBatchfileMalformed, self.options_parser._genomes_to_process, '', path_batchfile,
                              'fna')
        finally:
            shutil.rmtree(tmp_genome_dir)

    def test__genomes_to_process__batchfile__blank_genome_id(self):
        """ Test that a batchfile containing a blank genome id throws an exception. """
        try:
            tmp_genome_dir = tempfile.mkdtemp()
            path_batchfile = os.path.join(tmp_genome_dir, 'batchfile.txt')
            path_genome_1 = os.path.join(tmp_genome_dir, 'genome_1.fna')
            path_genome_2 = os.path.join(tmp_genome_dir, 'genome_2.fna')
            open(path_genome_1, 'a').close()
            open(path_genome_2, 'a').close()

            with open(path_batchfile, 'a') as f:
                f.write('%s\tgenome_1\n' % path_genome_1)
                f.write('\n')
                f.write('%s\t\n' % path_genome_2)

            self.assertRaises(GenomeBatchfileMalformed, self.options_parser._genomes_to_process, '', path_batchfile,
                              'fna')
        finally:
            shutil.rmtree(tmp_genome_dir)

    def test__genomes_to_process__batchfile__duplicate_genome_id(self):
        """ Test that a batchfile containing duplicate genome ids throws an exception. """
        # Branch 1: The number of columns are not equal to 2.
        try:
            tmp_genome_dir = tempfile.mkdtemp()
            path_batchfile = os.path.join(tmp_genome_dir, 'batchfile.txt')
            path_genome_1 = os.path.join(tmp_genome_dir, 'genome_1.fna')
            path_genome_2 = os.path.join(tmp_genome_dir, 'genome_2.fna')
            open(path_genome_1, 'a').close()
            open(path_genome_2, 'a').close()

            with open(path_batchfile, 'a') as f:
                f.write('%s\tgenome_1\n' % path_genome_1)
                f.write('\n')
                f.write('%s\tgenome_1\n' % path_genome_2)

            self.assertRaises(GenomeBatchfileMalformed, self.options_parser._genomes_to_process, '', path_batchfile,
                              'fna')
        finally:
            shutil.rmtree(tmp_genome_dir)

    def test__genomes_to_process__batchfile__invalid_genome_id(self):
        """ Test that a batchfile containing duplicate genome ids throws an exception. """
        # Branch 1: The number of columns are not equal to 2.
        try:
            tmp_genome_dir = tempfile.mkdtemp()
            path_batchfile_1 = os.path.join(tmp_genome_dir, 'batchfile_1.txt')
            path_batchfile_2 = os.path.join(tmp_genome_dir, 'batchfile_2.txt')
            path_batchfile_3 = os.path.join(tmp_genome_dir, 'batchfile_3.txt')
            path_genome_1 = os.path.join(tmp_genome_dir, 'genome_1.fna')
            path_genome_2 = os.path.join(tmp_genome_dir, 'genome_2.fna')
            open(path_genome_1, 'a').close()
            open(path_genome_2, 'a').close()

            with open(path_batchfile_1, 'a') as f:
                f.write('%s\tgenome_1\n' % path_genome_1)
                f.write('\n')
                f.write('%s\tGB_genome_2\n' % path_genome_2)

            with open(path_batchfile_2, 'a') as f:
                f.write('%s\tgenome_1\n' % path_genome_1)
                f.write('\n')
                f.write('%s\tRS_genome_2\n' % path_genome_2)

            with open(path_batchfile_3, 'a') as f:
                f.write('%s\tgenome_1\n' % path_genome_1)
                f.write('\n')
                f.write('%s\tUBAgenome_2\n' % path_genome_2)

            self.assertRaises(GenomeNameInvalid, self.options_parser._genomes_to_process, '', path_batchfile_1, 'fna')
            self.assertRaises(GenomeNameInvalid, self.options_parser._genomes_to_process, '', path_batchfile_2, 'fna')
            self.assertRaises(GenomeNameInvalid, self.options_parser._genomes_to_process, '', path_batchfile_3, 'fna')
        finally:
            shutil.rmtree(tmp_genome_dir)

    def test__genomes_to_process__no_files(self):
        """ Test that an exception is thrown if no files are found to process """
        # Branch 1 : genome_dir is specified
        try:
            tmp_genome_dir = tempfile.mkdtemp()
            self.assertRaises(NoGenomesFound, self.options_parser._genomes_to_process, tmp_genome_dir, '', 'fna')
        finally:
            shutil.rmtree(tmp_genome_dir)

        # Branch 2: batchfile is specified
        try:
            tmp_genome_dir = tempfile.mkdtemp()
            path_batchfile = os.path.join(tmp_genome_dir, 'batchfile.txt')
            open(path_batchfile, 'a').close()
            self.assertRaises(NoGenomesFound, self.options_parser._genomes_to_process, '', path_batchfile, 'fna')
        finally:
            shutil.rmtree(tmp_genome_dir)

    def test__marker_set_id(self):
        """ Test that the correct marker set id is returned """
        self.assertEqual(self.options_parser._marker_set_id(True, False, False), 'bac120')
        self.assertEqual(self.options_parser._marker_set_id(False, True, False), 'ar122')
        self.assertEqual(self.options_parser._marker_set_id(False, False, True), 'rps23')

    def test_identify__genome_dir_raises_io_exception(self):
        """ Test that the identify method raises an exception on invalid genome_dir """
        options = argparse.ArgumentParser()
        options.genome_dir = os.path.join(tempfile.gettempdir(), 'non-existent-dir')
        self.assertRaises(BioLibDirNotFound, self.options_parser.identify, options)

    def test_identify__batchfile_raises_io_exception(self):
        """ Test that the identify method raises an exception on invalid batchfile """
        options = argparse.ArgumentParser()
        options.genome_dir = None
        options.batchfile = os.path.join(tempfile.gettempdir(), 'non-existent-file.txt')
        self.assertRaises(BioLibFileNotFound, self.options_parser.identify, options)

    def test_align__identify_dir_raises_io_exception(self):
        """ Test that the align method raises an exception on invalid identify dir """
        options = argparse.ArgumentParser()
        options.identify_dir = os.path.join(tempfile.gettempdir(), 'non-existent-dir')
        self.assertRaises(BioLibDirNotFound, self.options_parser.align, options)

    def test_infer__msa_raises_io_exception(self):
        """ Test that the infer method raises an exception on invalid MSA """
        options = argparse.ArgumentParser()
        options.msa_file = os.path.join(tempfile.gettempdir(), 'non-existent-msa.txt')
        self.assertRaises(BioLibFileNotFound, self.options_parser.infer, options)

    def test_classify__align_dir_raises_io_exception(self):
        """ Test that the classify method raises an exception on invalid align dir """
        options = argparse.ArgumentParser()
        options.align_dir = os.path.join(tempfile.gettempdir(), 'non-existent-dir')
        self.assertRaises(BioLibDirNotFound, self.options_parser.classify, options)

    def test_root__no_tree_raises_io_exception(self):
        """ Test that the infer method raises an exception on invalid tree """
        options = argparse.ArgumentParser()
        options.input_tree = os.path.join(tempfile.gettempdir(), 'non-existent-tree.tree')
        self.assertRaises(BioLibFileNotFound, self.options_parser.root, options)

    def test_decorate__no_tree_raises_io_exception(self):
        """ Test that the infer method raises an exception on invalid tree """
        options = argparse.ArgumentParser()
        options.input_tree = os.path.join(tempfile.gettempdir(), 'non-existent-tree.tree')
        self.assertRaises(BioLibFileNotFound, self.options_parser.decorate, options)

    def test_trim_msa__mask_file(self):
        """ Test that the expected result is returned when running trim_msa with mask_file """
        try:
            tmp_out_dir = tempfile.mkdtemp()
            path_untrimmed_msa = os.path.join(tmp_out_dir, 'untrimmed_msa.fasta')
            path_mask_file = os.path.join(tmp_out_dir, 'mask_file.txt')
            path_output = os.path.join(tmp_out_dir, 'trimmed_msa.fasta')

            with open(path_untrimmed_msa, 'w') as f:
                f.write('>genome_1\n')
                f.write('ALGPVW\n')
                f.write('>genome_2\n')
                f.write('WVPGLA\n')

            with open(path_mask_file, 'w') as f:
                f.write('010010\n')

            options = argparse.ArgumentParser()
            # Required arguments
            options.untrimmed_msa = path_untrimmed_msa
            options.output = path_output
            # Mutex arguments
            options.mask_file = path_mask_file
            options.reference_mask = None

            self.options_parser.trim_msa(options)

            results = dict()
            with open(path_output, 'r') as f:
                re_hits = re.findall(r'>(.+)\n(.+)\n', f.read())
                for gid, seq in re_hits:
                    results[gid] = seq

            expected = {'genome_1': 'LV', 'genome_2': 'VL'}

            self.assertDictEqual(results, expected)

            pass
        finally:
            shutil.rmtree(tmp_out_dir)

    def test_trim_msa__reference_mask_arc(self):
        """ Test that the expected result is returned when running trim_msa with archaeal reference_mask """
        try:
            tmp_out_dir = tempfile.mkdtemp()
            path_untrimmed_msa = os.path.join(tmp_out_dir, 'untrimmed_msa.fasta')
            path_output = os.path.join(tmp_out_dir, 'trimmed_msa.fasta')

            msa_str = str()
            while len(msa_str) < 32675:
                msa_str += 'ALGPVW'
            msa_str = msa_str[0:32675]

            with open(path_untrimmed_msa, 'w') as f:
                f.write('>genome_1\n')
                f.write('%s\n' % msa_str)
                f.write('>genome_2\n')
                f.write('%s\n' % msa_str[::-1])

            options = argparse.ArgumentParser()
            # Required arguments
            options.untrimmed_msa = path_untrimmed_msa
            options.output = path_output
            # Mutex arguments
            options.mask_file = None
            options.reference_mask = 'arc'

            self.options_parser.trim_msa(options)

            results = dict()
            with open(path_output, 'r') as f:
                re_hits = re.findall(r'>(.+)\n(.+)\n', f.read())
                for gid, seq in re_hits:
                    results[gid] = hashlib.sha256(seq.encode()).hexdigest()

            expected = {'genome_1': '332b8cd125a36c375196064e136efab78db38e41bbd8bd8484243531bc57df6d',
                        'genome_2': '84e91b9f5fa1ec0bedc0097233044e6dd0e79557bb6df3625928dc9573795989'}

            self.assertDictEqual(results, expected)
        finally:
            shutil.rmtree(tmp_out_dir)

    def test_trim_msa__reference_mask_bac(self):
        """ Test that the expected result is returned when running trim_msa with bacterial reference_mask """
        try:
            tmp_out_dir = tempfile.mkdtemp()
            path_untrimmed_msa = os.path.join(tmp_out_dir, 'untrimmed_msa.fasta')
            path_output = os.path.join(tmp_out_dir, 'trimmed_msa.fasta')

            msa_str = str()
            while len(msa_str) < 41155:
                msa_str += 'ALGPVW'
            msa_str = msa_str[0:41155]

            with open(path_untrimmed_msa, 'w') as f:
                f.write('>genome_1\n')
                f.write('%s\n' % msa_str)
                f.write('>genome_2\n')
                f.write('%s\n' % msa_str[::-1])

            options = argparse.ArgumentParser()
            # Required arguments
            options.untrimmed_msa = path_untrimmed_msa
            options.output = path_output
            # Mutex arguments
            options.mask_file = None
            options.reference_mask = 'bac'

            self.options_parser.trim_msa(options)

            results = dict()
            with open(path_output, 'r') as f:
                re_hits = re.findall(r'>(.+)\n(.+)\n', f.read())
                for gid, seq in re_hits:
                    results[gid] = hashlib.sha256(seq.encode()).hexdigest()

            expected = {'genome_1': '35e080f9ab7d318e8f4a7cef46ce6044bd9c538e6fbe8a69b17431df44bd5a81',
                        'genome_2': 'bb4beed69063dad0092a809ee5854ff124da0b55c651edd50c47b1d8fdff0d7b'}

            self.assertDictEqual(results, expected)
        finally:
            shutil.rmtree(tmp_out_dir)

    def test_export_msa__arc(self):
        """ Test that the untrimmed archaeal MSA is exported correctly """
        try:
            tmp_out_dir = tempfile.mkdtemp()
            path_out = os.path.join(tmp_out_dir, 'output.fasta')

            options = argparse.ArgumentParser()
            options.domain = 'arc'
            options.output = path_out

            self.options_parser.export_msa(options)

            with open(path_out, 'rb') as f:
                out_hash = hashlib.sha256(f.read()).hexdigest()
            self.assertEqual(out_hash, '11eb12b91ab20c43824abafb909ccc20bed84a8609a9bf82748b2cdbdd8b7aad')
        finally:
            shutil.rmtree(tmp_out_dir)

    def test_export_msa__bac(self):
        """ Test that the untrimmed bacterial MSA is exported correctly """
        try:
            tmp_out_dir = tempfile.mkdtemp()
            path_out = os.path.join(tmp_out_dir, 'output.fasta')

            options = argparse.ArgumentParser()
            options.domain = 'bac'
            options.output = path_out

            self.options_parser.export_msa(options)

            with open(path_out, 'rb') as f:
                out_hash = hashlib.sha256(f.read()).hexdigest()
            self.assertEqual(out_hash, '50dde1e96df9533def7c7047a1e8627d4ad566db10f8ab3de72751e62c4ac10a')
        finally:
            shutil.rmtree(tmp_out_dir)

