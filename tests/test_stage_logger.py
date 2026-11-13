import json
import os
import shutil
import tempfile
import unittest

from gtdbtk.files.stage_logger import ClassifyStep, StageLogger


class TestStageLogger(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix='gtdbtk_tmp_')
        self.sl = StageLogger()
        self.sl.path = os.path.join(self.tmp, 'gtdbtk.json')
        self.sl.steps = []

    def tearDown(self):
        self.sl.steps = []
        shutil.rmtree(self.tmp)

    def test_classify_step_legacy_pplacer_keys(self):
        # gtdbtk.json written by GTDB-Tk < 2.9 (pplacer).
        legacy = {'steps': [{'name': 'classify', 'status': 'completed', 'align_dir': '/x/align',
                             'scratch_dir': '/x/scratch', 'full_tree': True}]}
        with open(self.sl.path, 'w') as fh:
            json.dump(legacy, fh)
        self.sl.read_existing_steps()
        step = self.sl.get_stage(ClassifyStep)
        self.assertTrue(step.is_complete())
        self.assertEqual(step.align_dir, '/x/align')

        self.sl.write()
        with open(self.sl.path) as fh:
            written = json.load(fh)['steps'][0]
        self.assertNotIn('scratch_dir', written)
        self.assertNotIn('full_tree', written)


if __name__ == '__main__':
    unittest.main()
