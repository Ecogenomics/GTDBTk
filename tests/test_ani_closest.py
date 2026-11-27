import shutil
import tempfile
import unittest
from unittest import mock

from gtdbtk import ani_rep as ani_rep_mod
from gtdbtk.ani_rep import ANIClosestFile
from gtdbtk.biolib_lite.common import canonical_gid


class _Radii:
    def __init__(self, radii=None):
        self.radii = radii or {}

    def get_rep_ani(self, gid):
        return self.radii.get(gid, 95.0)


class TestANIClosestFile(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix='gtdbtk_ani_closest_')

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def test_no_result_rows_match_header_column_count(self):
        ref = 'GCF_000000001.1'
        taxonomy = {
            canonical_gid(ref): [
                'd__Bacteria', 'p__P', 'c__C', 'o__O', 'f__F', 'g__G', 's__G sp1'
            ]
        }
        genomes = {
            'g_hit': '/tmp/g_hit.fna',
            'g_low_af': '/tmp/g_low_af.fna',
            'g_no_skani_hit': '/tmp/g_no_skani_hit.fna',
        }
        results = {
            'g_hit': {ref: {'ani': 96.5, 'af': 0.8}},
            'g_low_af': {ref: {'ani': 96.5, 'af': 0.1}},
        }

        with mock.patch.object(ani_rep_mod, 'GTDBRadiiFile', _Radii):
            closest_file = ANIClosestFile(
                self.tmp, 'pfx', results, genomes, min_af=0.5, taxonomy=taxonomy
            )

        with open(closest_file.path) as fh:
            lines = [line.rstrip('\n').split('\t') for line in fh]

        header = lines[0]
        self.assertEqual(len(header), 6)
        rows = {row[0]: row for row in lines[1:]}
        self.assertEqual(set(rows), set(genomes))

        for gid, row in rows.items():
            self.assertEqual(
                len(row),
                len(header),
                f'Row for {gid} has {len(row)} columns, expected {len(header)}: {row}',
            )

        self.assertEqual(rows['g_hit'][1], ref)
        self.assertEqual(rows['g_hit'][5], 'True')
        self.assertEqual(rows['g_low_af'][1:], ['no result'] * 5)
        self.assertEqual(rows['g_no_skani_hit'][1:], ['no result'] * 5)


if __name__ == '__main__':
    unittest.main()
