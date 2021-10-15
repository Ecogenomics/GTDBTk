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
import shutil
import tempfile
import unittest

from gtdbtk.biolib_lite.exceptions import *
from gtdbtk.exceptions import *
from gtdbtk.main import OptionsParser


class TestOptionsParser(unittest.TestCase):

    def setUp(self):
        self.options_parser = OptionsParser('-1')
        self.dir_tmp = tempfile.mkdtemp(prefix='gtdbtk_tmp_')
        pass

    def tearDown(self):
        shutil.rmtree(self.dir_tmp)

    def test__verify_genome_id__valid(self):
        """ Test that a valid genome id returns True. """
        self.assertTrue(self.options_parser._verify_genome_id('genome_1'))

    def test__verify_genome_id__invalid(self):
        """ Test that invalid genome ids throw an exception. """
        for c in list('()[],;='):
            self.assertRaises(GTDBTkExit, self.options_parser._verify_genome_id, 'genome%s1' % c)

    def test__genomes_to_process__genome_dir__valid(self):
        """ Test that the expected results are returned when using genome_dir. """
        open(os.path.join(self.dir_tmp, 'genome_1.fna'), 'a').close()
        open(os.path.join(self.dir_tmp, 'genome_2.fna'), 'a').close()
        open(os.path.join(self.dir_tmp, 'other_file.txt'), 'a').close()
        results, tln_table = self.options_parser._genomes_to_process(self.dir_tmp, '', 'fna')
        expected = {'genome_1': os.path.join(self.dir_tmp, 'genome_1.fna'),
                    'genome_2': os.path.join(self.dir_tmp, 'genome_2.fna')}
        self.assertDictEqual(results, expected)

    def test__genomes_to_process__batchfile__valid(self):
        """ Test that the expected results are returned when using batchfile """
        path_batchfile = os.path.join(self.dir_tmp, 'batchfile.txt')
        path_genome_1 = os.path.join(self.dir_tmp, 'genome_1.fna')
        path_genome_2 = os.path.join(self.dir_tmp, 'genome_2.fna')
        open(path_genome_1, 'a').close()
        open(path_genome_2, 'a').close()

        with open(path_batchfile, 'a') as f:
            f.write(f'{path_genome_1}\tgenome_1\n')
            f.write('\n')
            f.write(f'{path_genome_2}\tgenome_2\t4\n')

        results, tln_table = self.options_parser._genomes_to_process('', path_batchfile, 'fna')
        expected = {'genome_1': path_genome_1, 'genome_2': path_genome_2}
        expected_tln = {'genome_2': 4}
        self.assertDictEqual(results, expected)
        self.assertDictEqual(tln_table, expected_tln)

    def test__genomes_to_process__batchfile__invalid_columns(self):
        """ Test that a batchfile containing columns not equal to 2 throws an exception. """
        path_batchfile = os.path.join(self.dir_tmp, 'batchfile.txt')
        path_genome_1 = os.path.join(self.dir_tmp, 'genome_1.fna')
        path_genome_2 = os.path.join(self.dir_tmp, 'genome_2.fna')
        open(path_genome_1, 'a').close()
        open(path_genome_2, 'a').close()

        with open(path_batchfile, 'a') as f:
            f.write('%s\tgenome_1\n' % path_genome_1)
            f.write('\n')
            f.write('%s\tgenome_2\tfoo\n' % path_genome_2)

        self.assertRaises(GTDBTkExit, self.options_parser._genomes_to_process, '', path_batchfile,
                          'fna')

    def test__genomes_to_process__batchfile__blank_genome_path(self):
        """ Test that a batchfile containing a blank genome path throws an exception. """
        path_batchfile = os.path.join(self.dir_tmp, 'batchfile.txt')
        path_genome_1 = os.path.join(self.dir_tmp, 'genome_1.fna')
        path_genome_2 = os.path.join(self.dir_tmp, 'genome_2.fna')
        open(path_genome_1, 'a').close()
        open(path_genome_2, 'a').close()

        with open(path_batchfile, 'a') as f:
            f.write('%s\tgenome_1\n' % path_genome_1)
            f.write('\n')
            f.write('%s\tgenome_2\n' % '')

        self.assertRaises(GTDBTkExit, self.options_parser._genomes_to_process, '', path_batchfile,
                          'fna')

    def test__genomes_to_process__batchfile__blank_genome_id(self):
        """ Test that a batchfile containing a blank genome id throws an exception. """
        path_batchfile = os.path.join(self.dir_tmp, 'batchfile.txt')
        path_genome_1 = os.path.join(self.dir_tmp, 'genome_1.fna')
        path_genome_2 = os.path.join(self.dir_tmp, 'genome_2.fna')
        open(path_genome_1, 'a').close()
        open(path_genome_2, 'a').close()

        with open(path_batchfile, 'a') as f:
            f.write('%s\tgenome_1\n' % path_genome_1)
            f.write('\n')
            f.write('%s\t\n' % path_genome_2)

        self.assertRaises(GTDBTkExit, self.options_parser._genomes_to_process, '', path_batchfile,
                          'fna')

    def test__genomes_to_process__batchfile__duplicate_genome_id(self):
        """ Test that a batchfile containing duplicate genome ids throws an exception. """
        # Branch 1: The number of columns are not equal to 2.
        path_batchfile = os.path.join(self.dir_tmp, 'batchfile.txt')
        path_genome_1 = os.path.join(self.dir_tmp, 'genome_1.fna')
        path_genome_2 = os.path.join(self.dir_tmp, 'genome_2.fna')
        open(path_genome_1, 'a').close()
        open(path_genome_2, 'a').close()

        with open(path_batchfile, 'a') as f:
            f.write('%s\tgenome_1\n' % path_genome_1)
            f.write('\n')
            f.write('%s\tgenome_1\n' % path_genome_2)

        self.assertRaises(GTDBTkExit, self.options_parser._genomes_to_process, '', path_batchfile,
                          'fna')

    # def test__genomes_to_process__batchfile__invalid_genome_id(self):
    #     """ Test that a batchfile containing duplicate genome ids throws an exception. """
    #     # Branch 1: The number of columns are not equal to 2.
    #     path_batchfile_1 = os.path.join(self.dir_tmp, 'batchfile_1.txt')
    #     path_batchfile_2 = os.path.join(self.dir_tmp, 'batchfile_2.txt')
    #     path_batchfile_3 = os.path.join(self.dir_tmp, 'batchfile_3.txt')
    #     path_genome_1 = os.path.join(self.dir_tmp, 'genome_1.fna')
    #     path_genome_2 = os.path.join(self.dir_tmp, 'genome_2.fna')
    #     open(path_genome_1, 'a').close()
    #     open(path_genome_2, 'a').close()
    #
    #     with open(path_batchfile_1, 'a') as f:
    #         f.write('%s\tgenome_1\n' % path_genome_1)
    #         f.write('\n')
    #         f.write('%s\tGB_genome_2\n' % path_genome_2)
    #
    #     with open(path_batchfile_2, 'a') as f:
    #         f.write('%s\tgenome_1\n' % path_genome_1)
    #         f.write('\n')
    #         f.write('%s\tRS_genome_2\n' % path_genome_2)
    #
    #     with open(path_batchfile_3, 'a') as f:
    #         f.write('%s\tgenome_1\n' % path_genome_1)
    #         f.write('\n')
    #         f.write('%s\tUBAgenome_2\n' % path_genome_2)
    #
    #     self.assertRaises(GTDBTkExit, self.options_parser._genomes_to_process, '', path_batchfile_1, 'fna')
    #     self.assertRaises(GTDBTkExit, self.options_parser._genomes_to_process, '', path_batchfile_2, 'fna')
    #     self.assertRaises(GTDBTkExit, self.options_parser._genomes_to_process, '', path_batchfile_3, 'fna')

    def test__genomes_to_process__no_files(self):
        """ Test that an exception is thrown if no files are found to process """
        # Branch 1 : genome_dir is specified
        tmp_genome_dir = tempfile.mkdtemp()
        try:
            self.assertRaises(GTDBTkExit, self.options_parser._genomes_to_process, tmp_genome_dir, '', 'fna')
        finally:
            shutil.rmtree(tmp_genome_dir)

        # Branch 2: batchfile is specified
        tmp_genome_dir = tempfile.mkdtemp()
        try:
            path_batchfile = os.path.join(tmp_genome_dir, 'batchfile.txt')
            open(path_batchfile, 'a').close()
            self.assertRaises(GTDBTkExit, self.options_parser._genomes_to_process, '', path_batchfile, 'fna')
        finally:
            shutil.rmtree(tmp_genome_dir)

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

    def test_run_test(self):
        """Test that the user-test method runs correctly"""
        options = argparse.ArgumentParser()
        options.out_dir = self.dir_tmp
        options.cpus = 3
        self.assertTrue(self.options_parser.run_test(options))

    # def test_run_test__throws_exception(self):
    #     """Test that the user-test method fails correctly"""
    #     options = argparse.ArgumentParser()
    #     options.out_dir = self.dir_tmp
    #     os.mkdir(os.path.join(self.dir_tmp, 'genomes'))
    #     options.cpus = 3
    #     self.assertRaises(GTDBTkTestFailure, self.options_parser.run_test, options)

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
