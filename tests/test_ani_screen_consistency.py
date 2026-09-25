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

"""Regression tests for Issue #717: `classify` and `classify_wf` must give the same answer.

`classify` runs skani itself (Classify._sort_skani_results_pre_pplacer), while
`classify_wf` runs the ANI screen first (ANIScreener.sort_skani_ani_screen), writes
its results to disk and reloads them in classify. These tests check that both paths
assign the same species, list the same other related references and report the same
AF, and that the hits of unassigned genomes survive the round trip to disk.

No reference data or third-party binaries are needed: reference-data lookups are mocked.
"""

import copy
import json
import logging
import os
import random
import shutil
import tempfile
import unittest
from unittest import mock

from gtdbtk import ani_screen as ani_screen_mod
from gtdbtk.ani_rep import ANISummaryFile, SkaniRawHitsFile
from gtdbtk.biolib_lite.common import canonical_gid
from gtdbtk.classify import Classify
from gtdbtk.config.output import PATH_ANISCREEN_UNASSIGNED_HITS
from gtdbtk.external.skani import SkANI
from gtdbtk.files.stage_logger import ANIScreenStep, StageLogger
from gtdbtk.tools import add_ncbi_prefix


class _Radii(object):
    """Stand-in for GTDBRadiiFile: canonical gid -> radius (default 95)."""

    def __init__(self, radii=None):
        self.radii = radii or {}

    def get_rep_ani(self, gid):
        return self.radii.get(gid, 95.0)


def _taxonomy(refs, domain='d__Archaea'):
    """Canonical-keyed taxonomy (as used by the ANI screen)."""
    return {canonical_gid(r): [domain, 'p__', 'c__', 'o__', 'f__', 'g__', f's__{r}'] for r in refs}


def _other_ref_ids(other_refs):
    return [e.split(',')[0].strip() for e in other_refs.split(';')] if other_refs else []


class TestANIScreenConsistency(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix='gtdbtk_tmp_')
        self.af_threshold = 0.5

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    # ------------------------------------------------------------------ helpers
    def _screen(self, skani_results, taxonomy, radii=None):
        """classify_wf path: ANIScreener.sort_skani_ani_screen -> {gid: (ref, other_refs)}."""
        with mock.patch.object(ani_screen_mod, 'GTDBRadiiFile', lambda: _Radii(radii)), \
                mock.patch.object(ani_screen_mod.Classify, 'parse_radius_file', return_value={}):
            screener = ani_screen_mod.ANIScreener(1, self.af_threshold)
            out = screener.sort_skani_ani_screen(copy.deepcopy(skani_results), taxonomy)
        return {gid: (ref, hit.get('other_related_refs'))
                for dom in out.values() for gid, d in dom.items() for ref, hit in d.items()}

    def _classify(self, skani_results, taxonomy, radii=None):
        """classify path: Classify._sort_skani_results_pre_pplacer -> {gid: (ref, other_refs)}."""
        c = Classify.__new__(Classify)
        c.logger = logging.getLogger('timestamp')
        c.af_threshold = self.af_threshold
        c.gtdb_radii = _Radii(radii)
        c.species_radius = {}
        c.gtdb_taxonomy = {add_ncbi_prefix(r): v for r, v in
                           ((r, taxonomy[canonical_gid(r)]) for g in skani_results.values() for r in g)}
        out = c._sort_skani_results_pre_pplacer(copy.deepcopy(skani_results), {})
        return {row.gid: (row.closest_genome_ref, row.other_related_refs)
                for rows in out.values() for row in rows}

    # -------------------------------------------------------- species rule / refs
    def test_other_refs_af_filtered_and_identical(self):
        """AF < min_af hits are not listed; assigned rep is not repeated; both paths agree."""
        hits = {'GCA_018263495.1': {'ani': 96.12, 'af': 0.682},  # assigned
                'GCA_014384145.1': {'ani': 95.43, 'af': 0.68},
                'GCA_013390745.1': {'ani': 94.57, 'af': 0.738},
                'GCA_000402075.1': {'ani': 94.56, 'af': 0.773},
                'GCA_013390375.1': {'ani': 89.91, 'af': 0.557},
                'GCA_013390765.1': {'ani': 88.99, 'af': 0.309},  # AF < 0.5
                'GCA_902606965.1': {'ani': 87.64, 'af': 0.238}}  # AF < 0.5
        tax = _taxonomy(hits)
        wf, cl = self._screen({'g': hits}, tax), self._classify({'g': hits}, tax)
        self.assertEqual(wf, cl)
        ref, other = wf['g']
        self.assertEqual(ref, 'GCA_018263495.1')
        self.assertEqual(_other_ref_ids(other),
                         ['GCA_014384145.1', 'GCA_013390745.1', 'GCA_000402075.1', 'GCA_013390375.1'])

    def test_single_other_ref_kept(self):
        hits = {'GCA_000000001.1': {'ani': 97.0, 'af': 0.8},
                'GCA_000000002.1': {'ani': 94.0, 'af': 0.7},
                'GCA_000000003.1': {'ani': 93.0, 'af': 0.2}}
        tax = _taxonomy(hits)
        wf, cl = self._screen({'g': hits}, tax), self._classify({'g': hits}, tax)
        self.assertEqual(wf, cl)
        self.assertEqual(_other_ref_ids(wf['g'][1]), ['GCA_000000002.1'])

    def test_raised_radius_not_rescued_by_lower_hit(self):
        """Closest AF-passing hit fails its raised radius -> unassigned, even if a lower hit passes its own."""
        hits = {'GCA_000000001.1': {'ani': 97.0, 'af': 0.8},  # radius 98 -> fails
                'GCA_000000002.1': {'ani': 96.0, 'af': 0.8}}  # radius 95 -> would pass
        radii = {canonical_gid('GCA_000000001.1'): 98.0}
        tax = _taxonomy(hits, 'd__Bacteria')
        self.assertEqual(self._screen({'g': hits}, tax, radii), {})
        self.assertEqual(self._classify({'g': hits}, tax, radii), {})

    def test_other_refs_are_the_50_closest(self):
        hits = {f'GCF_{i:09d}.1': {'ani': 80.0 + i * 0.25, 'af': 0.9} for i in range(1, 61)}
        items = list(hits.items())
        random.Random(1).shuffle(items)
        hits = dict(items)
        tax = _taxonomy(hits, 'd__Bacteria')
        wf, cl = self._screen({'g': hits}, tax), self._classify({'g': hits}, tax)
        self.assertEqual(wf, cl)
        self.assertEqual(wf['g'][0], 'GCF_000000060.1')
        self.assertEqual(_other_ref_ids(wf['g'][1]), [f'GCF_{i:09d}.1' for i in range(59, 9, -1)])

    def test_randomised_equivalence(self):
        rnd = random.Random(7)
        refs = [f'GCA_{i:09d}.1' for i in range(1, 201)]
        radii = {canonical_gid(r): rnd.choice([95.0, 95.0, 95.0, 96.5, 97.5, 98.0]) for r in refs}
        tax = {canonical_gid(r): [rnd.choice(['d__Bacteria', 'd__Archaea']), 'p__', 'c__', 'o__', 'f__', 'g__',
                                  f's__{r}'] for r in refs}
        data = {f'u{g}': {r: {'ani': round(rnd.uniform(78, 100), 2), 'af': round(rnd.uniform(0.05, 1), 3)}
                          for r in rnd.sample(refs, rnd.randint(1, 60))} for g in range(2000)}
        wf, cl = self._screen(data, tax, radii), self._classify(data, tax, radii)
        self.assertGreater(len(wf), 0)
        self.assertEqual(wf, cl)

    # ------------------------------------------------------------------ AF rounding
    def test_af_rounding_identical_in_both_paths(self):
        """classify rounds AF once; classify_wf goes through ani_summary.tsv first. Same 3-digit AF."""
        rnd = random.Random(3)
        afs = [i / 100 for i in range(3000, 10001)] + [round(rnd.uniform(30, 100), 4) for _ in range(5000)]
        afs += [62.95, 62.9502, 89.45, 83.35]
        ref = 'GCA_000000009.1'
        sk = SkANI.__new__(SkANI)
        parsed = sk.parse_results(iter([(f'q{i}', ref, 97.0, a, 0.0) for i, a in enumerate(afs)]))
        classify_af = {q: round(h[ref]['af'], 3) for q, h in parsed.items()}

        f = ANISummaryFile(self.tmp, 'e2e', copy.deepcopy(parsed), _taxonomy([ref]), 'ar53')
        f.write(ani_screen_step=True)
        wf_af = {q: round(v[ref]['af'], 3) for q, v in ANISummaryFile(f.path).read().items()}
        self.assertEqual(classify_af, wf_af)
        self.assertEqual(parsed['q%d' % afs.index(89.45)][ref]['af'], 0.8945)

    # -------------------------------------------------- unassigned hits round trip
    def test_raw_hits_round_trip_exact(self):
        raw = {'u1': {'GCF_000000001.1': {'ani': 97.123456789, 'af': 0.81234567891234},
                      'GCA_000000002.2': {'ani': 94.99999999999999, 'af': 0.1 + 0.2}},
               'u2': {'GCF_000000003.1': {'ani': 80.0, 'af': 0.0}}}
        path = os.path.join(self.tmp, 'hits.tsv.gz')
        SkaniRawHitsFile.write(path, raw, {'u1', 'u2', 'not_in_results'})
        self.assertEqual(SkaniRawHitsFile.read(path), raw)
        SkaniRawHitsFile.write(path, raw, set())
        self.assertEqual(SkaniRawHitsFile.read(path), {})

    def test_run_aniscreen_writes_only_unassigned_genomes(self):
        raw = {'u1': {'GCF_000000001.1': {'ani': 97.1, 'af': 0.81}},    # assigned
               'u2': {'GCF_000000003.1': {'ani': 80.0, 'af': 0.6}},     # below radius
               'u3': {'GCF_000000004.1': {'ani': 99.0, 'af': 0.3}}}     # low AF
        tax = _taxonomy(['GCF_000000001.1', 'GCF_000000003.1', 'GCF_000000004.1'], 'd__Bacteria')
        with mock.patch.object(ani_screen_mod.ANIRep, 'run_skani', return_value=copy.deepcopy(raw)), \
                mock.patch.object(ani_screen_mod.Taxonomy, 'read', return_value=tax), \
                mock.patch.object(ani_screen_mod, 'CONFIG', mock.Mock(TAXONOMY_FILE='taxonomy.tsv')), \
                mock.patch.object(ani_screen_mod, 'GTDBRadiiFile', _Radii), \
                mock.patch.object(ani_screen_mod.Classify, 'parse_radius_file', return_value={}):
            classified, _reports, hits_file = ani_screen_mod.ANIScreener(1, self.af_threshold).run_aniscreen(
                {'u1': 'a', 'u2': 'b', 'u3': 'c'}, self.tmp, 'pfx')
        self.assertEqual(set(classified['bac120']), {'u1'})
        self.assertEqual(hits_file, os.path.join(self.tmp, PATH_ANISCREEN_UNASSIGNED_HITS.format(prefix='pfx')))
        self.assertEqual(SkaniRawHitsFile.read(hits_file), {'u2': raw['u2'], 'u3': raw['u3']})

    def test_stage_logger_field_and_legacy_json(self):
        sl = StageLogger()
        sl.path = os.path.join(self.tmp, 'gtdbtk.json')
        sl.steps = []
        step = ANIScreenStep()
        step.status = 'completed'
        step.unassigned_hits_file = '/x/hits.tsv.gz'
        sl.steps.append(step)
        sl.write()
        sl.steps = []
        sl.read_existing_steps()
        self.assertEqual(sl.get_stage(ANIScreenStep).unassigned_hits_file, '/x/hits.tsv.gz')

        with open(sl.path) as fh:  # gtdbtk.json written by GTDB-Tk < 2.8.0
            data = json.load(fh)
        del data['steps'][0]['unassigned_hits_file']
        with open(sl.path, 'w') as fh:
            json.dump(data, fh)
        sl.steps = []
        sl.read_existing_steps()
        self.assertIsNone(sl.get_stage(ANIScreenStep).unassigned_hits_file)

    def test_classify_run_reloads_hits_only_when_skani_skipped(self):
        class _Stop(Exception):
            pass

        def run(genes, path):
            c = Classify.__new__(Classify)
            c.logger = logging.getLogger('timestamp')
            c.cpus = 1
            seen = {}

            def fake_read(p):
                seen['read'] = p
                return {}

            with mock.patch.object(SkaniRawHitsFile, 'read', side_effect=fake_read), \
                    mock.patch('gtdbtk.classify.ClassifySummaryFileAR53', side_effect=_Stop):
                try:
                    c.run(genomes={}, align_dir=self.tmp, out_dir=self.tmp, prefix='pfx', skip_ani_screen=True,
                          genes=genes, ani_summary_files={}, unassigned_hits_file=path, all_classified_ani=True)
                except _Stop:
                    pass
            return seen

        hits_file = os.path.join(self.tmp, 'hits.tsv.gz')
        SkaniRawHitsFile.write(hits_file, {}, set())
        self.assertEqual(run(False, hits_file), {'read': hits_file})
        self.assertEqual(run(True, hits_file), {})       # --genes: no skani hits
        self.assertEqual(run(False, None), {})           # no ANI screen step
        self.assertEqual(run(False, os.path.join(self.tmp, 'missing.gz')), {})  # warns, does not crash


if __name__ == '__main__':
    unittest.main()
