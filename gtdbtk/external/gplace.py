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

import logging
import multiprocessing as mp
import os
import queue
import re
import subprocess
import sys

from tqdm import tqdm

from gtdbtk.exceptions import GplaceException, TogException
from gtdbtk.tools import get_proc_memory_gb


class Gplace(object):
    """Phylogenetic placement of genomes into a reference tree
    (https://github.com/cmc-aau/gPlace).
    """

    def __init__(self):
        """ Instantiate the class. """
        self.logger = logging.getLogger('timestamp')
        self.version = self._get_version()

    def _get_version(self):
        try:
            env = os.environ.copy()
            proc = subprocess.Popen(['gplace', '--version'], stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE, env=env, encoding='utf-8')

            output, error = proc.communicate()
            return output.strip()[1]
        except:
            return "(version unavailable)"

    def run(self, cpus, model, ref_pkg, json_out, msa_file, gplace_out,
            gplace_max_memory=None):
        """Place genomes into a reference tree.

        Args:
            cpus (int): The number of threads to use.
            model (str): The model to use. PROT: LG, WAG, JTT. NT: GTR.
            ref_pkg (str): The path to the reference package.
            json_out (str): The path to write the json output to.
            msa_file (str): The path to the input MSA file.
            gplace_out (str): Where to write the gplace output file.
            gplace_max_memory (int): The maximum memory to use in GB.

        Raises:
            GplaceException: if a non-zero exit code, or if the json output
                              file isn't generated.

        """
        args = ['gplace','place' '-m', model, '-j', str(cpus), '-c', ref_pkg, '-o',
                json_out, msa_file]
        if gplace_max_memory:
            args.append('--max-memory')
            args.append(str(gplace_max_memory))
        self.logger.debug(' '.join(args))

        out_q = mp.Queue()
        pid = mp.Value('i', 0)
        # print(f'Running gplace with the following command: {" ".join(args)}')
        # sys.exit(0)
        p_worker = mp.Process(target=self._worker, args=(
            args, out_q, gplace_out, pid))
        p_writer = mp.Process(target=self._writer, args=(out_q, pid))

        try:
            p_worker.start()
            p_writer.start()

            p_worker.join()
            out_q.put(None)
            p_writer.join()

            if p_worker.exitcode != 0:
                raise GplaceException(
                    'An error was encountered while running gplace.')
        except Exception:
            p_worker.terminate()
            p_writer.terminate()
            raise

        if not os.path.isfile(json_out):
            self.logger.error('gplace returned a zero exit code but no output '
                              'file was generated.')
            raise GplaceException

    def _worker(self, args, out_q, gplace_out, pid):
        """The worker thread writes the piped output of gplace to disk and
        shares it with the writer thread for logging."""
        with subprocess.Popen(args, stdout=subprocess.PIPE, encoding='utf-8') as proc:
            with pid.get_lock():
                pid.value = proc.pid

            with open(gplace_out, 'w') as fh:
                while True:
                    line = proc.stdout.readline()
                    if not line:
                        break
                    fh.write(f'{line}')
                    out_q.put(line)
            proc.wait()

            if proc.returncode != 0:
                raise GplaceException('An error was encountered while '
                                       'running gplace, check the log '
                                       'file: {}'.format(gplace_out))

    def tog(self, gplace_json_out, tree_file):
        """ Convert the gplace json output into a newick tree.

        Args:
            pplacer_json_out (str): The path to the output of pplacer.
            tree_file (str): The path to output the newick file to.

        Raises:
            GplaceException: If a non-zero exit code is returned, or the tree file
                          isn't output.
        """

        args = ['gplace', 'tog', '-o', tree_file, gplace_json_out]
        proc = subprocess.Popen(args, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
        proc_out, proc_err = proc.communicate()

        if proc.returncode != 0:
            self.logger.error('An error was encountered while running tog.')
            raise GplaceException(proc_err)

        if not os.path.isfile(gplace_json_out):
            self.logger.error('gplace tog returned a zero exit code but no output '
                              'file was generated.')
            raise GplaceException
