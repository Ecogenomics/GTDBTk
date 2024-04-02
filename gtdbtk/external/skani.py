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
import re
import shutil
import subprocess
import tempfile

from gtdbtk.exceptions import GTDBTkExit
from gtdbtk.tools import tqdm_log


class SkANI(object):
    """Python wrapper for skani (https://github.com/bluenote-1577/skani)"""

    def __init__(self, cpus, force_single):
        """Instantiate the class.

        Parameters
        ----------
        cpus : int
            The number of CPUs to use.
        force_single : bool
            True if the ql and rl calls should be done individually.
        """
        self.cpus = max(cpus, 1)
        self.force_single = force_single
        self.logger = logging.getLogger('timestamp')
        self.version = self._get_version()
        self._suppress_v1_warning = False

    @staticmethod
    def _get_version():
        """Returns the version of skani on the system path.

        Returns
        -------
        str
            The string containing the skani version.
        """
        try:
            proc = subprocess.Popen(['skani', '-V'], stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE, encoding='utf-8')
            stdout, stderr = proc.communicate()
            version = re.search(r'skani (.+)', stdout)
            return version.group(1)
        except Exception as e:
            return 'unknown'


    def run(self, dict_compare, dict_paths,preset='--medium', report_progress=True):
        """Runs skani in batch mode.

        Parameters
        ----------
        dict_compare : dict[str, set[str]]
            All query to reference comparisons to be made.
        dict_paths : dict[str, str]
            The path for each genome id being compared.

        Returns
        -------
        dict[str, dict[str, dict[str, float]]]
            A dictionary containing the ANI and AF for each comparison."""

        # Create the multiprocessing items.
        manager = mp.Manager()
        q_worker = manager.Queue()
        q_writer = manager.Queue()
        q_results = manager.Queue()

        # Populate the queue of comparisons in forwards and reverse direction.
        n_total = 0
        list_comparisons = [(user_id, ref_id) for user_id, values in dict_compare.items() for ref_id in values]
        for qry, ref in list_comparisons:
            n_total += 1
            q_worker.put((qry, ref, dict_paths.get(qry), dict_paths.get(ref), preset))

        # Set the terminate condition for each worker thread.
        [q_worker.put(None) for _ in range(self.cpus)]

        # Create each of the processes
        p_workers = [mp.Process(target=self._worker,
                                args=(q_worker, q_writer, q_results))
                     for _ in range(self.cpus)]

        p_writer = mp.Process(target=self._writer, args=(q_writer, n_total))

        # Start each of the threads.
        try:
            # Start the writer and each processing thread.
            p_writer.start()
            for p_worker in p_workers:
                p_worker.start()

            # Wait until each worker has finished.
            for p_worker in p_workers:
                p_worker.join()

                # Gracefully terminate the program.
                if p_worker.exitcode != 0:
                    raise GTDBTkExit('skani returned a non-zero exit code.')

            # Stop the writer thread.
            q_writer.put(None)
            p_writer.join()

        except Exception:
            for p in p_workers:
                p.terminate()
            p_writer.terminate()
            raise

        # Process and return each of the results obtained
        q_results.put(None)
        return self._parse_result_queue(q_results)

    def _worker(self, q_worker, q_writer, q_results):
        """Operates skani in list mode.

        Parameters
        ----------
        q_worker : mp.Queue
            A multiprocessing queue containing the available jobs.
        q_writer : mp.Queue
            A multiprocessing queue to track progress.
        q_results : mp.Queue
            A multiprocessing queue containing raw results.
        """
        while True:
            # Retrieve the next item, stop if the sentinel is found.
            job = q_worker.get(block=True)
            if job is None:
                break

            # Extract the values
            q, r, q_path,r_path,preset = job

            # Run skani
            result = self.run_proc(q, r, q_path, r_path, preset)
            q_results.put(result)
            q_writer.put(True)

        return True

    def _writer(self, q_writer, n_total):
        """The writer function, which reports the progress of the workers.

        Parameters
        ----------
        q_writer : mp.Queue
            A queue of genome ids which have been processed.
        n_total : int
            The total number of items to be processed.
        """
        with tqdm_log(unit='comparison', total=n_total) as p_bar:
            for _ in iter(q_writer.get, None):
                p_bar.update()

    def run_proc(self, qid, rid, ql, rl, preset):
        """Runs the skani process.

        Parameters
        ----------
        q : str
            The path to the query genome.
        r : str
            The path to the reference genome.
        ql : str
            The path to the query list file.
        rl : str
            The path to the reference list file.
        output : str
            The path to the output file.

        Returns
        -------
        dict[str, dict[str, float]]
            The ANI/AF of the query genomes to the reference genomes.
        """
        args = ['skani', 'dist',
                preset,
                '-q', ql,
                '-r', rl,
                '-o', '/dev/stdout']
        # self.logger.debug(' '.join(args))
        proc = subprocess.Popen(args, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE, encoding='utf-8')
        stdout, stderr = proc.communicate()

        if proc.returncode != 0:
            self.logger.error('STDOUT:\n' + stdout)
            self.logger.error('STDERR:\n' + stderr)
            raise GTDBTkExit('skani returned a non-zero exit code.')

        result_lines = stdout.splitlines()
        if len(result_lines) == 2:
            tokens = result_lines[1].split('\t')
            ani = float(tokens[2])
            af_r = float(tokens[3])
            af_q = float(tokens[4])
            ani_af = (qid, rid, ani, af_r, af_q)
        else:
        #     # genomes too divergent to determine ANI and AF
        #     # with skani so default to zeros
            ani_af = "null"
        #     ani_af = (qid, rid, 0.0, 0.0, 0.0)


        return ani_af


    def _maybe_write_list(self, d_genomes, path):
        """Writes a query/reference list to disk.

        Parameters
        ----------
        d_genomes : dict[str, str]
            A dictionary containing the key as the accession and value as path.
        path : str
            The path to write the file to.
        """
        if d_genomes is None or path is None:
            return
        with open(path, 'w') as fh:
            for gid, gid_path in d_genomes.items():
                fh.write(f'{gid_path}\n')

    def _parse_result_queue(self, q_results):
        """Creates the output dictionary given the results from skani

        Parameters
        ----------
        q_results : mp.Queue
            A multiprocessing queue containing raw results.

        Returns
        -------
        dict[str, dict[str, dict[str, float]]]
            The ANI/AF of the query genome to all reference genomes.
        """
        out = dict()
        while True:
            job = q_results.get(block=True)
            if job == 'null':
                continue
            if job is None:
                break

            qid, rid, ani, af_r, af_q = job
            max_af = max(af_r, af_q)
            # af is a percent, we need to divide it by 100
            max_af = max_af/100
            if qid not in out:
                out[qid] = {rid: {'ani': ani, 'af': max_af}}
            else:
                out[qid][rid] = {'ani': ani, 'af': max_af}

        return out
