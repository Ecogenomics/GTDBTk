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

import multiprocessing as mp
import os
import subprocess
import sys
import unittest

# Runs Parallel in a fresh interpreter, since the start method can only be set
# once per process. If the manager process dies on start, the parent waits on
# the consumer queue forever, so a timeout is how the test sees the failure.
SCRIPT = """
import multiprocessing as mp
import sys

mp.set_start_method(sys.argv[1])

from gtdbtk.biolib_lite.parallel import Parallel

total = Parallel(cpus=2).run(abs, lambda x, acc: (acc or 0) + x, [-1, -2, 3])
print(total)
"""


class TestBioLibParallel(unittest.TestCase):

    def test_run_all_start_methods(self):
        """ Parallel.run works under fork, spawn and forkserver """
        env = dict(os.environ, PYTHONPATH=os.pathsep.join(sys.path))
        for method in mp.get_all_start_methods():
            with self.subTest(method=method):
                proc = subprocess.run([sys.executable, '-c', SCRIPT, method],
                                      env=env, stdout=subprocess.PIPE,
                                      stderr=subprocess.PIPE,
                                      universal_newlines=True, timeout=60)
                self.assertEqual(proc.returncode, 0, proc.stderr)
                self.assertEqual(proc.stdout.strip(), '6')
