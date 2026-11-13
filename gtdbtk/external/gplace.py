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
import threading

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
        args = ['gplace','place', '-m', model, '-j', str(cpus), '-c', ref_pkg, '-o',
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
        # Thread, so gplace warnings reach the GTDB-Tk logger under spawn (macOS).
        p_writer = threading.Thread(target=self._writer, args=(out_q, pid), daemon=True)

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
            out_q.put(None)
            raise

        if not os.path.isfile(json_out):
            self.logger.error('gplace returned a zero exit code but no output '
                              'file was generated.')
            raise GplaceException

    def _worker(self, args, out_q, gplace_out, pid):
        """The worker thread writes the piped output of gplace to disk and
        shares it with the writer thread for logging."""
        with subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                              encoding='utf-8') as proc:
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

    def _writer(self, out_q, pid):
        """The writer thread reports gplace progress from the piped log lines,
        and memory usage while the reference tree is cached."""
        states = ['Starting gplace',
                  'Reading user alignment',
                  'Reading reference package',
                  'Caching likelihood vectors',
                  'Ranking edges',
                  'Placing genomes',
                  'Writing placements']
        cur_state, peak_gb = 0, None
        bar_fmt = '==> {desc}.'
        with tqdm(bar_format=bar_fmt, miniters=1, leave=False) as p_bar:
            p_bar.set_description_str(desc=f'Step 1 of {len(states)}: {states[0]}')
            while True:
                try:
                    line = out_q.get(block=True, timeout=5)
                    if not line:
                        break
                    # gplace lines: [YYYY-MM-DD HH:MM:SS] LEVEL: message
                    hit = re.match(r'\[[^\]]+\] (\w+): (.*)', line.strip())
                    if not hit:
                        continue
                    level, msg = hit.groups()
                    if level in ('WARNING', 'ERROR'):
                        self.logger.log(logging.getLevelName(level), f'gplace: {msg}')
                        continue

                    if msg.startswith('Read '):
                        cur_state = 2
                    elif re.match(r'Placing [\d,]+ queries on', msg):
                        cur_state = 3
                        peak = re.search(r'estimated peak memory ([\d.]+) GiB', msg)
                        peak_gb = float(peak.group(1)) if peak else None
                    elif msg.startswith(('Cached likelihood vectors', 'Boundary vectors')):
                        cur_state = 4
                    elif msg.startswith(('Scored ', 'Ranked ')):
                        cur_state = 5
                    elif msg.startswith(('Refined the kept placements', 'Wrote ')):
                        cur_state = 6
                    elif cur_state == 0:
                        cur_state = 1

                    desc = f'Step {cur_state + 1} of {len(states)}: {states[cur_state]}'
                    # "Prefilter: tile 3/12", "Placed 800/1,000 queries.", "Baseball round 2: 800/1,000 ..."
                    count = re.search(r'([\d,]+)/([\d,]+)', msg)
                    if count and cur_state in (4, 5):
                        n_done, n_total = (int(x.replace(',', '')) for x in count.groups())
                        unit = 'tile' if cur_state == 4 else 'genome'
                        desc += f': {unit} {n_done:,} of {n_total:,} ({n_done / n_total:.2%})'
                    p_bar.set_description_str(desc=desc)

                # Caching is the long silent step: report memory against gplace's estimate.
                except queue.Empty:
                    if cur_state == 3:
                        virt, res = get_proc_memory_gb(pid.value)
                        if res:
                            est = f' of ~{peak_gb:.1f} GiB estimated' if peak_gb else ''
                            p_bar.set_description_str(desc=f'Step 4 of {len(states)}: {states[3]} '
                                                           f'({res:.2f} GB{est})')
                except Exception:
                    pass

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
