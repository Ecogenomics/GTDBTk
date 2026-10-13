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

"""Hit ordering: `>`, `max()` and `sorted()` must agree (Hit.__gt__ used to `raise NotImplemented`)."""

import random
import unittest

from gtdbtk.files.marker.tophit import Hit


class TestHitOrdering(unittest.TestCase):

    def test_comparison_operators(self):
        weak = Hit('g1', 'TIGR00001', 1e-5, 50.0)
        strong = Hit('g2', 'TIGR00002', 1e-20, 120.0)
        self.assertTrue(strong > weak)
        self.assertTrue(weak < strong)
        self.assertTrue(strong >= weak and weak <= strong)
        self.assertFalse(weak > strong)
        self.assertEqual(max([weak, strong]), strong)

    def test_tie_breaks(self):
        # same bit-score: lower e-value wins; then lower hmm id; then lower gene id
        self.assertGreater(Hit('g', 'H1', 1e-9, 80.0), Hit('g', 'H1', 1e-5, 80.0))
        self.assertGreater(Hit('g', 'H1', 1e-9, 80.0), Hit('g', 'H2', 1e-9, 80.0))
        self.assertGreater(Hit('a', 'H1', 1e-9, 80.0), Hit('b', 'H1', 1e-9, 80.0))

    def test_max_matches_sorted(self):
        """max() (now used for top hits) picks the same hit as the previous sorted(reverse=True)[0]."""
        rnd = random.Random(11)
        for _ in range(2000):
            hits = [Hit(rnd.choice('abc'), rnd.choice(['H1', 'H2', 'H3']),
                        rnd.choice([1e-30, 1e-10, 1e-5]), rnd.choice([40.0, 80.0, 120.0]))
                    for _ in range(rnd.randint(1, 8))]
            self.assertEqual(max(hits), sorted(hits, reverse=True)[0])


if __name__ == '__main__':
    unittest.main()
