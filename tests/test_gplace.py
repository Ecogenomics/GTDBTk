import io
import os
import shutil
import tempfile
import unittest
from unittest import mock

from gtdbtk.biolib_lite.logger import logger_setup
from gtdbtk.cli import get_main_parser
from gtdbtk.exceptions import GplaceException, GTDBTkExit
from gtdbtk.external.gplace import Gplace
from gtdbtk.main import OptionsParser


class _FakePopen(object):
    """Stands in for gplace: records each command and writes its -o file."""
    calls = []
    version_out = 'gPlace v0.1.1\n'
    write_output = True

    def __init__(self, args, **_kwargs):
        _FakePopen.calls.append(args)
        self.returncode, self.pid = 0, 1
        if self.write_output and '-o' in args:
            open(args[args.index('-o') + 1], 'w').close()
        self.stdout = io.StringIO('[2026-11-13 10:00:00] INFO: Read 2 sequences from user.msa\n')

    def communicate(self):
        return self.version_out, ''

    def wait(self):
        return self.returncode

    def __enter__(self):
        return self

    def __exit__(self, *_args):
        return False


class _SyncProcess(object):
    """mp.Process run in the test process, so the Popen mock applies to the worker."""

    def __init__(self, target, args):
        self.target, self.args, self.exitcode = target, args, None

    def start(self):
        try:
            self.target(*self.args)
            self.exitcode = 0
        except Exception:
            self.exitcode = 1

    def join(self):
        pass

    def terminate(self):
        pass


@mock.patch('gtdbtk.external.gplace.mp.Process', _SyncProcess)
@mock.patch('gtdbtk.external.gplace.subprocess.Popen', _FakePopen)
class TestGplace(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix='gtdbtk_tmp_')
        _FakePopen.calls, _FakePopen.write_output = [], True
        _FakePopen.version_out = 'gPlace v0.1.1\n'

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def test_version(self):
        self.assertEqual(Gplace().version, '0.1.1')
        _FakePopen.version_out = 'gplace: version 0.1.0 Copyright 2026 Pierre-Alain Chaumeil\n'
        self.assertEqual(Gplace().version, '0.1.0')

    def test_run_command(self):
        json_out = os.path.join(self.tmp, 'gplace.json')
        gplace_out = os.path.join(self.tmp, 'gplace.out')
        Gplace().run(4, 'wag', 'ref.refpkg', json_out, 'user.msa', gplace_out, '64G')
        self.assertEqual(_FakePopen.calls[-1],
                         ['gplace', 'place', '-m', 'wag', '-j', '4', '-c', 'ref.refpkg', '-o', json_out,
                          'user.msa', '--max_memory', '64G'])
        with open(gplace_out) as fh:
            self.assertIn('Read 2 sequences', fh.read())

    def test_run_without_max_memory(self):
        json_out = os.path.join(self.tmp, 'gplace.json')
        Gplace().run(1, 'wag', 'ref.refpkg', json_out, 'user.msa', os.path.join(self.tmp, 'gplace.out'))
        self.assertNotIn('--max_memory', _FakePopen.calls[-1])

    def test_run_no_output(self):
        _FakePopen.write_output = False
        with self.assertRaises(GplaceException):
            Gplace().run(1, 'wag', 'ref.refpkg', os.path.join(self.tmp, 'gplace.json'), 'user.msa',
                         os.path.join(self.tmp, 'gplace.out'))

    def test_tog(self):
        tree_file = os.path.join(self.tmp, 'classify.tree')
        Gplace().tog('gplace.json', tree_file)
        self.assertEqual(_FakePopen.calls[-1], ['gplace', 'tog', '-o', tree_file, 'gplace.json'])

        _FakePopen.write_output = False
        with self.assertRaises(GplaceException):
            Gplace().tog('gplace.json', os.path.join(self.tmp, 'missing.tree'))


@mock.patch.object(OptionsParser, '_check_package_compatibility', lambda self: None)
@mock.patch.object(OptionsParser, 'classify')
class TestGplaceOptions(unittest.TestCase):

    def setUp(self):
        logger_setup(None, "gtdbtk.log", "GTDB-Tk", 'unittest', True)
        self.parser = get_main_parser()
        self.base = ['classify', '--genome_dir', 'g', '--align_dir', 'a', '--out_dir', 'o']

    def test_deprecated_pplacer_flags(self, classify):
        options = self.parser.parse_args(self.base + ['-f', '--scratch_dir', 's', '--pplacer_cpus', '4',
                                                      '--skip_pplacer'])
        with self.assertLogs('timestamp', level='WARNING') as logs:
            OptionsParser('unittest').parse_options(options)
        self.assertTrue(classify.called)
        self.assertEqual(options.gplace_cpus, 4)
        self.assertTrue(options.skip_gplace)
        self.assertEqual(sum('ignored with gplace' in x for x in logs.output), 2)

    def test_gplace_cpus_wins_over_pplacer_cpus(self, _classify):
        options = self.parser.parse_args(self.base + ['--pplacer_cpus', '4', '--gplace_cpus', '8'])
        OptionsParser('unittest').parse_options(options)
        self.assertEqual(options.gplace_cpus, 8)

    def test_gplace_max_memory(self, classify):
        for value in ('64G', '500M', '1.5T', '64', '64GiB'):
            options = self.parser.parse_args(self.base + ['--gplace_max_memory', value])
            OptionsParser('unittest').parse_options(options)
        self.assertEqual(classify.call_count, 5)
        options = self.parser.parse_args(self.base + ['--gplace_max_memory', '64X'])
        with self.assertRaises(GTDBTkExit):
            OptionsParser('unittest').parse_options(options)


if __name__ == '__main__':
    unittest.main()
