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

from tqdm import tqdm

from gtdbtk.exceptions import PplacerException, TogException
from gtdbtk.tools import get_proc_memory_gb


class Pplacer(object):
    """Phylogenetic placement of genomes into a reference tree
    (http://matsen.fredhutch.org/pplacer/).
    """

    def __init__(self):
        """ Instantiate the class. """
        self.logger = logging.getLogger('timestamp')
        self.version = self._get_version()

    def _get_version(self):
        try:
            env = os.environ.copy()
            proc = subprocess.Popen(['pplacer', '--version'], stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE, env=env, encoding='utf-8')

            output, error = proc.communicate()
            return output.strip()
        except:
            return "(version unavailable)"

    def run(self, cpus, model, ref_pkg, json_out, msa_file, pplacer_out,
            mmap_file=None):
        """Place genomes into a reference tree.

        Args:
            cpus (int): The number of threads to use.
            model (str): The model to use. PROT: LG, WAG, JTT. NT: GTR.
            ref_pkg (str): The path to the reference package.
            json_out (str): The path to write the json output to.
            msa_file (str): The path to the input MSA file.
            pplacer_out (str): Where to write the pplacer output file.
            mmap_file (str, optional): The path to write a scratch file to.

        Raises:
            PplacerException: if a non-zero exit code, or if the json output
                              file isn't generated.

        """
        args = ['pplacer', '-m', model, '-j', str(cpus), '-c', ref_pkg, '-o',
                json_out, msa_file]
        if mmap_file:
            args.append('--mmap-file')
            args.append(mmap_file)
        self.logger.debug(' '.join(args))

        out_q = mp.Queue()
        pid = mp.Value('i', 0)
        p_worker = mp.Process(target=self._worker, args=(
            args, out_q, pplacer_out, pid))
        p_writer = mp.Process(target=self._writer, args=(out_q, pid))

        try:
            p_worker.start()
            p_writer.start()

            p_worker.join()
            out_q.put(None)
            p_writer.join()

            if p_worker.exitcode != 0:
                raise PplacerException(
                    'An error was encountered while running pplacer.')
        except Exception:
            p_worker.terminate()
            p_writer.terminate()
            raise
        finally:
            if mmap_file:
                os.remove(mmap_file)

        if not os.path.isfile(json_out):
            self.logger.error('pplacer returned a zero exit code but no output '
                              'file was generated.')
            raise PplacerException

    def _worker(self, args, out_q, pplacer_out, pid):
        """The worker thread writes the piped output of pplacer to disk and
        shares it with the writer thread for logging."""
        with subprocess.Popen(args, stdout=subprocess.PIPE, encoding='utf-8') as proc:
            with pid.get_lock():
                pid.value = proc.pid

            with open(pplacer_out, 'w') as fh:
                while True:
                    line = proc.stdout.readline()
                    if not line:
                        break
                    fh.write(f'{line}')
                    out_q.put(line)
            proc.wait()

            if proc.returncode != 0:
                raise PplacerException('An error was encountered while '
                                       'running pplacer, check the log '
                                       'file: {}'.format(pplacer_out))

    def _writer(self, out_q, pid):
        """The writer subprocess is able to report on newly piped events from
        subprocess in the worker thread, and report on memory usage while
        waiting for new commands."""
        states = ['Reading user alignment',
                  'Reading reference alignment',
                  'Pre-masking sequences',
                  'Determining figs',
                  'Allocating memory for internal nodes',
                  'Caching likelihood information on reference tree',
                  'Pulling exponents',
                  'Preparing the edges for baseball',
                  'Placing genomes']
        cur_state = None
        n_total, n_placed = None, 0
        # Create a progress bar which tracks the initialization step of
        # pplacer.
        bar_fmt = '==> {desc}.'
        with tqdm(bar_format=bar_fmt, miniters=1, leave=False) as p_bar:
            p_bar.set_description_str(desc='Step 1 of 9: Starting pplacer')
            while True:
                try:
                    state = out_q.get(block=True, timeout=5)
                    if not state:
                        break
                    elif state.startswith('Running pplacer'):
                        cur_state = 0
                    elif state.startswith("Didn't find any reference"):
                        cur_state = 1
                    elif state.startswith('Pre-masking sequences'):
                        cur_state = 2
                    elif state.startswith('Determining figs'):
                        cur_state = 3
                    elif state.startswith('Allocating memory for internal'):
                        cur_state = 4
                    elif state.startswith('Caching likelihood information'):
                        cur_state = 5
                    elif state.startswith('Pulling exponents'):
                        cur_state = 6
                    elif state.startswith('Preparing the edges'):
                        cur_state = 7
                    elif state.startswith('working on '):
                        cur_state = 8
                    else:
                        cur_state = None

                    # Error, display the log file line.
                    if not cur_state:
                        p_bar.set_description_str(desc=state.strip())

                    # New line written for each placed genome, update the
                    # count.
                    elif cur_state == 8:
                        if not n_total:
                            n_total = int(
                                re.search(r'\((\d+)\/(\d+)\)', state).group(2))
                        n_placed += 1
                        p_bar.set_description_str(desc=f'Step 9 of {len(states)}: placing '
                                                  f'genome {n_placed:,} of {n_total:,} '
                                                  f'({n_placed / n_total:.2%})')

                    # Display the output for the current state.
                    else:
                        p_bar.set_description_str(desc=f'Step {cur_state + 1} of {len(states)}: '
                                                  f'{states[cur_state + 1]}')

                # Report the memory usage if at a memory-reportable state.
                except queue.Empty:
                    if cur_state == 3:
                        virt, res = get_proc_memory_gb(pid.value)
                        if virt and res:
                            p_bar.set_description_str(desc=f'Step {cur_state + 1} of {len(states)} '
                                                      f'{states[4]} ({virt:.2f} GB)')
                        else:
                            p_bar.set_description_str(desc=f'Step {cur_state + 1} of {len(states)} '
                                                           f'{states[4]}')
                    elif cur_state == 4:
                        virt, res = get_proc_memory_gb(pid.value)
                        if virt and res:
                            p_bar.set_description_str(desc=f'Step {cur_state + 1} of {len(states)}: '
                                                      f'{states[5]} ({res:.2f}/{virt:.2f} GB, '
                                                      f'{res / virt:.2%})')
                        else:
                            p_bar.set_description_str(desc=f'Step {cur_state + 1} of {len(states)}: '
                                                           f'{states[5]}')
                except Exception:
                    pass

    def tog(self, pplacer_json_out, tree_file):
        """ Convert the pplacer json output into a newick tree.

        Args:
            pplacer_json_out (str): The path to the output of pplacer.
            tree_file (str): The path to output the newick file to.

        Raises:
            TogException: If a non-zero exit code is returned, or the tree file
                          isn't output.
        """

        args = ['guppy', 'tog', '-o', tree_file, pplacer_json_out]
        proc = subprocess.Popen(args, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
        proc_out, proc_err = proc.communicate()

        if proc.returncode != 0:
            self.logger.error('An error was encountered while running tog.')
            raise TogException(proc_err)

        if not os.path.isfile(pplacer_json_out):
            self.logger.error('tog returned a zero exit code but no output '
                              'file was generated.')
            raise TogException
