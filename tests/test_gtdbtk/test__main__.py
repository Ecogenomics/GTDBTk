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
import subprocess
import tempfile
import unittest

from gtdbtk.config.config import CONCAT_AR53, CONCAT_BAC120, MASK_AR53, MASK_DIR, MASK_BAC120
from gtdbtk.files.classify_summary import ClassifySummaryFileAR53
from gtdbtk.tools import sha256


class TestMain(unittest.TestCase):
    """Used to test the CLI. Currently just checks all exit codes."""

    genome_dir = 'gtdbtk/tests/data/genomes/'
    genome_dir_gz = 'tests/data/genomes_gz/'
    cpus = '30'

    def setUp(self):
        """Create a new temporary directory for each method."""
        self.dir_tmp_obj = tempfile.TemporaryDirectory(prefix='gtdbtk_tmp')
        self.dir_tmp = self.dir_tmp_obj.name

    def tearDown(self):
        """Cleanup after each method finishes."""
        self.dir_tmp_obj.cleanup()

    def test_classify_wf(self):
        args = ['python', '-m', 'gtdbtk', 'classify_wf', '--genome_dir',
                self.genome_dir, '--out_dir', self.dir_tmp, '--cpus', self.cpus]
        p = subprocess.Popen(args)
        p.wait()
        self.assertEqual(p.returncode, 0)

        classify_file = ClassifySummaryFileAR53(out_dir=self.dir_tmp, prefix='gtdbtk')
        classify_file.read()

        expected = {
            'genome_1': 'd__Archaea;p__Methanobacteriota;c__Methanobacteria;o__Methanobacteriales;f__Methanobacteriaceae;g__Methanobrevibacter;s__Methanobrevibacter ruminantium',
            'genome_2': 'd__Archaea;p__Thermoplasmatota;c__Thermoplasmata;o__Methanomassiliicoccales;f__Methanomethylophilaceae;g__VadinCA11;s__VadinCA11 sp002498365',
            'genome_3': 'd__Archaea;p__Thermoplasmatota;c__Thermoplasmata;o__Methanomassiliicoccales;f__Methanomethylophilaceae;g__VadinCA11;s__VadinCA11 sp002498365'
        }
        actual = {k: v.classification for k, v in classify_file.rows.items()}
        self.assertDictEqual(expected, actual)


    def test_de_novo_wf(self):
        args = ['python', '-m', 'gtdbtk', 'de_novo_wf', '--archaea', '--genome_dir',
                self.genome_dir, '--out_dir', self.dir_tmp, '--cpus', self.cpus,
                '--outgroup_taxon', 'p__Altarchaeota']
        proc = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='utf-8')
        stdout, stderr = proc.communicate()
        self.assertEqual(proc.returncode, 0)

    def test_identify(self):
        """Duplicated in classify"""
        pass

    def test_identify_gzipped_genomes(self):
        """Gene calling results should match if genomes are gzipped."""
        args = ['python', '-m', 'gtdbtk', 'identify', '--genome_dir',
                self.genome_dir_gz, '--out_dir', self.dir_tmp, '--cpus', self.cpus,
                '--extension', 'gz']
        p = subprocess.Popen(args)
        p.wait()
        self.assertEqual(p.returncode, 0)

    def test_align(self):
        """Duplicated in classify"""
        pass

    def test_classify(self):
        identify_dir = os.path.join(self.dir_tmp, 'identify')
        args = ['python', '-m', 'gtdbtk', 'identify', '--genome_dir',
                self.genome_dir, '--out_dir', identify_dir, '--cpus', self.cpus]
        proc = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='utf-8')
        stdout, stderr = proc.communicate()
        self.assertEqual(proc.returncode, 0)

        align_dir = os.path.join(self.dir_tmp, 'align')
        args = ['python', '-m', 'gtdbtk', 'align', '--identify_dir',
                identify_dir, '--out_dir', align_dir, '--cpus', self.cpus]
        proc = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='utf-8')
        stdout, stderr = proc.communicate()
        self.assertEqual(proc.returncode, 0)

        classify_dir = os.path.join(self.dir_tmp, 'classify')
        args = ['python', '-m', 'gtdbtk', 'classify', '--genome_dir',
                self.genome_dir, '--align_dir', align_dir, '--out_dir', classify_dir,
                '--cpus', self.cpus]
        proc = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='utf-8')
        stdout, stderr = proc.communicate()
        self.assertEqual(proc.returncode, 0)

    def test_infer(self):
        """Duplicated in decorate"""
        pass

    def test_root(self):
        """Duplicated in decorate"""
        pass

    def test_decorate(self):
        identify_dir = os.path.join(self.dir_tmp, 'identify')
        args = ['python', '-m', 'gtdbtk', 'identify', '--genome_dir',
                self.genome_dir, '--out_dir', identify_dir, '--cpus', self.cpus]
        proc = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='utf-8')
        stdout, stderr = proc.communicate()
        self.assertEqual(proc.returncode, 0)

        align_dir = os.path.join(self.dir_tmp, 'align')
        args = ['python', '-m', 'gtdbtk', 'align', '--identify_dir',
                identify_dir, '--out_dir', align_dir, '--cpus', self.cpus]
        proc = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='utf-8')
        stdout, stderr = proc.communicate()
        self.assertEqual(proc.returncode, 0)

        msa_file = os.path.join(align_dir, 'align', 'gtdbtk.ar53.msa.fasta')
        infer_dir = os.path.join(self.dir_tmp, 'infer')
        args = ['python', '-m', 'gtdbtk', 'infer', '--msa_file',
                msa_file, '--out_dir', infer_dir, '--cpus', self.cpus]
        proc = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='utf-8')
        stdout, stderr = proc.communicate()
        self.assertEqual(proc.returncode, 0)

        input_tree = os.path.join(infer_dir, 'infer', 'intermediate_results', 'gtdbtk.unrooted.tree')
        rooted_tree = os.path.join(self.dir_tmp, 'root.tree')
        args = ['python', '-m', 'gtdbtk', 'root', '--input_tree',
                input_tree, '--outgroup_taxon', 'p__Altarchaeota',
                '--output_tree', rooted_tree]
        proc = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='utf-8')
        stdout, stderr = proc.communicate()
        self.assertEqual(proc.returncode, 0)

        decorated_tree = os.path.join(self.dir_tmp, 'decorated.tree')
        args = ['python', '-m', 'gtdbtk', 'decorate', '--input_tree',
                rooted_tree, '--output_tree', decorated_tree]
        proc = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='utf-8')
        stdout, stderr = proc.communicate()
        self.assertEqual(proc.returncode, 0)

        # Test infer ranks
        path_infer_tree = os.path.join(self.dir_tmp, 'infer.tree')
        args = ['python', '-m', 'gtdbtk', 'infer_ranks', '--input_tree', decorated_tree,
                '--ingroup_taxon', 'p__Altarchaeota', '--output_tree', path_infer_tree]
        proc = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='utf-8')
        stdout, stderr = proc.communicate()
        self.assertEqual(proc.returncode, 0)

    def test_infer_ranks(self):
        """Duplicated in decorate."""
        pass

    def test_ani_rep(self):
        args = ['python', '-m', 'gtdbtk', 'ani_rep', '--genome_dir', self.genome_dir,
                '--out_dir', self.dir_tmp, '--cpus', self.cpus]
        proc = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='utf-8')
        stdout, stderr = proc.communicate()
        self.assertEqual(proc.returncode, 0)

    def test_trim_msa(self):
        path_out = os.path.join(self.dir_tmp, 'msa.faa')
        args = ['python', '-m', 'gtdbtk', 'trim_msa', '--untrimmed_msa', CONCAT_AR53,
                '--output', path_out, '--reference_mask', 'arc']
        proc = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='utf-8')
        stdout, stderr = proc.communicate()
        self.assertEqual(proc.returncode, 0)

    def test_trim_msa_reference_mask_arc(self):
        with open(os.path.join(MASK_DIR, MASK_AR53), 'r') as f:
            mask = [x == '1' for x in f.read().strip()]

        path_untrimmed_msa = os.path.join(self.dir_tmp, 'untrimmed_msa.fasta')
        path_output = os.path.join(self.dir_tmp, 'trimmed_msa.fasta')

        choices = 'ARNDCQEGHILKMFPSTWYV-'
        expected = list()
        with open(path_untrimmed_msa, 'w') as f:
            f.write('>test123\n')
            for cur_char in mask:
                choice = random.choice(choices)
                if cur_char:
                    expected.append(choice)
                f.write(choice)

        args = ['python', '-m', 'gtdbtk', 'trim_msa', '--untrimmed_msa', path_untrimmed_msa,
                '--reference_mask', 'arc', '--output', path_output]
        p = subprocess.Popen(args)
        p.wait()
        self.assertEqual(p.returncode, 0)

        results = dict()
        with open(path_output, 'r') as f:
            import re
            re_hits = re.findall(r'>(.+)\n(.+)\n', f.read())
            for gid, seq in re_hits:
                results[gid] = seq

        expected = {'test123': ''.join(expected)}

        self.assertDictEqual(results, expected)

    def test_trim_msa_reference_mask_bac(self):
        with open(os.path.join(MASK_DIR, MASK_BAC120), 'r') as f:
            mask = [x == '1' for x in f.read().strip()]

        path_untrimmed_msa = os.path.join(self.dir_tmp, 'untrimmed_msa.fasta')
        path_output = os.path.join(self.dir_tmp, 'trimmed_msa.fasta')

        choices = 'ARNDCQEGHILKMFPSTWYV-'
        expected = list()
        with open(path_untrimmed_msa, 'w') as f:
            f.write('>test123\n')
            for cur_char in mask:
                choice = random.choice(choices)
                if cur_char:
                    expected.append(choice)
                f.write(choice)

        args = ['python', '-m', 'gtdbtk', 'trim_msa', '--untrimmed_msa', path_untrimmed_msa,
                '--reference_mask', 'bac', '--output', path_output]
        p = subprocess.Popen(args)
        p.wait()
        self.assertEqual(p.returncode, 0)

        results = dict()
        with open(path_output, 'r') as f:
            import re
            re_hits = re.findall(r'>(.+)\n(.+)\n', f.read())
            for gid, seq in re_hits:
                results[gid] = seq

        expected = {'test123': ''.join(expected)}

        self.assertDictEqual(results, expected)

    def test_trim_msa_mask_file(self):
        path_untrimmed_msa = os.path.join(self.dir_tmp, 'untrimmed_msa.fasta')
        path_mask_file = os.path.join(self.dir_tmp, 'mask_file.txt')
        path_output = os.path.join(self.dir_tmp, 'trimmed_msa.fasta')

        with open(path_untrimmed_msa, 'w') as f:
            f.write('>genome_1\n')
            f.write('ALGPVW\n')
            f.write('>genome_2\n')
            f.write('WVPGLA\n')

        with open(path_mask_file, 'w') as f:
            f.write('010010\n')

        args = ['python', '-m', 'gtdbtk', 'trim_msa', '--untrimmed_msa', path_untrimmed_msa,
                '--mask_file', path_mask_file, '--output', path_output]
        p = subprocess.Popen(args)
        p.wait()
        self.assertEqual(p.returncode, 0)

        results = dict()
        with open(path_output, 'r') as f:
            import re
            re_hits = re.findall(r'>(.+)\n(.+)\n', f.read())
            for gid, seq in re_hits:
                results[gid] = seq

        expected = {'genome_1': 'LV', 'genome_2': 'VL'}

        self.assertDictEqual(results, expected)

    def test_export_msa_arc(self):
        """Test that the MSA can be exported when using the CLI."""
        path_output = os.path.join(self.dir_tmp, 'msa.faa')
        args = ['python', '-m', 'gtdbtk', 'export_msa', '--domain', 'arc', '--output', path_output]
        p = subprocess.Popen(args)
        p.wait()
        self.assertEqual(p.returncode, 0)

        test_hash = sha256(path_output)
        true_hash = sha256(CONCAT_AR53)
        self.assertEqual(test_hash, true_hash)

    def test_export_msa_bac(self):
        """Test that the MSA can be exported when using the CLI."""
        path_output = os.path.join(self.dir_tmp, 'msa.faa')
        args = ['python', '-m', 'gtdbtk', 'export_msa', '--domain', 'bac', '--output', path_output]
        p = subprocess.Popen(args)
        p.wait()
        self.assertEqual(p.returncode, 0)

        test_hash = sha256(path_output)
        true_hash = sha256(CONCAT_BAC120)
        self.assertEqual(test_hash, true_hash)

    def test_test(self):
        args = ['python', '-m', 'gtdbtk', 'test', '--out_dir', self.dir_tmp, '--cpus', self.cpus]
        proc = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='utf-8')
        stdout, stderr = proc.communicate()
        self.assertEqual(proc.returncode, 0)

    def test_check_install(self):
        args = ['python', '-m', 'gtdbtk', 'check_install']
        proc = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='utf-8')
        stdout, stderr = proc.communicate()
        self.assertEqual(proc.returncode, 0)
