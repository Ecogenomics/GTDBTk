import multiprocessing as mp
import os
import queue
import shutil
import tempfile
import types
import unittest
from unittest import mock

import gtdbtk.external.prodigal as prodigal_mod
from gtdbtk.external.prodigal import Prodigal


class TestProdigalForce(unittest.TestCase):
    """With --force, genomes Prodigal cannot process are reported as failed, with the same
    tuple shape as any other genome."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix='gtdbtk_tmp_')
        with mock.patch.object(Prodigal, '_get_version', lambda self: 'test'):
            self.prodigal = Prodigal(1, os.path.join(self.tmp, 'failed.tsv'), self.tmp,
                                     '_protein.faa', '_protein.fna', '_protein.gff', True)

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def _run(self, biolib_result):
        os.makedirs(os.path.join(self.tmp, 'g1'), exist_ok=True)
        fake = mock.Mock()
        fake.return_value.run.return_value = biolib_result
        with mock.patch.object(prodigal_mod, 'BioLibProdigal', fake):
            return self.prodigal._run_prodigal('g1', 'g1.fna', None)

    def test_header_only_fasta(self):
        # BioLib Prodigal returns nothing on an error (e.g. '>empty' without sequence)
        result = self._run({})
        self.assertEqual(len(result), 8)
        self.assertEqual(result[5:], (False, True, 'Prodigal failed to process the genome'))

    def test_empty_fasta(self):
        stats = types.SimpleNamespace(is_empty=True)
        result = self._run({'g1': stats})
        self.assertEqual(len(result), 8)
        self.assertEqual(result[5:], (False, True, 'Empty file'))      # not counted as skipped

    def test_worker_records_failed_genome(self):
        out_dict, worker_queue, writer_queue = {}, queue.Queue(), queue.Queue()
        n_skipped = mp.Value('i', 0)
        worker_queue.put(('g1', 'g1.fna', None))
        worker_queue.put(None)
        fake = mock.Mock()
        fake.return_value.run.return_value = {}
        with mock.patch.object(prodigal_mod, 'BioLibProdigal', fake):
            self.prodigal._worker(out_dict, worker_queue, writer_queue, n_skipped)
        self.assertEqual(n_skipped.value, 0)
        self.assertTrue(out_dict['g1']['is_empty'])
        self.assertEqual(out_dict['g1']['fail_reason'], 'Prodigal failed to process the genome')


if __name__ == '__main__':
    unittest.main()
