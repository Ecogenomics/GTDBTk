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
import shutil
import subprocess
import sys
import tempfile

from gtdbtk.config.config import MP_TIMEOUT
from gtdbtk.exceptions import GTDBTkExit


class FastANI(object):
    """Python wrapper for FastANI (https://github.com/ParBLiSS/FastANI)"""

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

    def run(self, dict_compare, dict_paths):
        """Runs FastANI in batch mode.

        Parameters
        ----------
        dict_compare : Dict[str, Dict[str, str]]
            All query to reference comparisons to be made.
        dict_paths : Dict[str, str]
            The path for each genome id being compared.

        Returns
        -------
        Dict[str, Dict[str, Dict[str, float]]]
            A dictionary containing the ANI and AF for each comparison."""

        # Create the multiprocessing items.
        manager = mp.Manager()
        q_worker = manager.Queue()
        q_writer = manager.Queue()
        q_results = manager.Queue()

        # Populate the queue of comparisons in forwards and reverse direction.
        n_total = 0
        if self.force_single:
            for qry_gid, ref_set in dict_compare.items():
                qry_path = dict_paths[qry_gid]

                for ref_gid in ref_set:
                    ref_path = dict_paths[ref_gid]

                    fwd_dict = {'q': dict(), 'r': dict(), 'qry': qry_gid}
                    rev_dict = {'q': dict(), 'r': dict(), 'qry': qry_gid}

                    fwd_dict['q'][qry_gid] = qry_path
                    fwd_dict['r'][ref_gid] = ref_path

                    rev_dict['q'][ref_gid] = ref_path
                    rev_dict['r'][qry_gid] = qry_path

                    q_worker.put(fwd_dict, timeout=MP_TIMEOUT)
                    q_worker.put(rev_dict, timeout=MP_TIMEOUT)
                    n_total += 2
        else:
            for qry_gid, ref_set in dict_compare.items():
                fwd_dict = {'ql': dict(), 'rl': dict(), 'qry': qry_gid}
                rev_dict = {'ql': dict(), 'rl': dict(), 'qry': qry_gid}

                qry_path = dict_paths[qry_gid]
                fwd_dict['ql'][qry_gid] = qry_path
                rev_dict['rl'][qry_gid] = qry_path

                for ref_gid in ref_set:
                    ref_path = dict_paths[ref_gid]
                    fwd_dict['rl'][ref_gid] = ref_path
                    rev_dict['ql'][ref_gid] = ref_path

                q_worker.put(fwd_dict, timeout=MP_TIMEOUT)
                q_worker.put(rev_dict, timeout=MP_TIMEOUT)
                n_total += 2

        # Set the terminate condition for each worker thread.
        [q_worker.put(None, timeout=MP_TIMEOUT) for _ in range(self.cpus)]

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
                p_worker.join(timeout=MP_TIMEOUT)

                # Gracefully terminate the program.
                if p_worker.exitcode != 0:
                    raise GTDBTkExit('FastANI returned a non-zero exit code.')

            # Stop the writer thread.
            q_writer.put(None, timeout=MP_TIMEOUT)
            p_writer.join(timeout=MP_TIMEOUT)

        except Exception:
            for p in p_workers:
                p.terminate()
            p_writer.terminate()
            raise

        # Process and return each of the results obtained
        path_to_gid = {v: k for k, v in dict_paths.items()}
        q_results.put(None, timeout=MP_TIMEOUT)
        return self._parse_result_queue(q_results, path_to_gid)

    def _worker(self, q_worker, q_writer, q_results):
        """Operates FastANI in list mode.

        Parameters
        ----------
        q_worker : Queue
            A multiprocessing queue containing the available jobs.
        q_writer : Queue
            A multiprocessing queue to track progress.
        q_results : Queue
            A multiprocessing queue containing raw results.
        """
        while True:
            # Retrieve the next item, stop if the sentinel is found.
            job = q_worker.get(block=True, timeout=MP_TIMEOUT)
            if job is None:
                break

            # Extract the values
            q = list(job.get('q').values())[0] if job.get('q') is not None else None
            r = list(job.get('r').values())[0] if job.get('r') is not None else None
            ql = job.get('ql')
            rl = job.get('rl')

            # Create a temporary directory to write the lists to.
            dir_tmp = tempfile.mkdtemp(prefix='gtdbtk_fastani_tmp_')
            path_qry = os.path.join(dir_tmp, 'ql.txt') if ql is not None else None
            path_ref = os.path.join(dir_tmp, 'rl.txt') if rl is not None else None
            path_out = os.path.join(dir_tmp, 'output.txt')

            try:
                # Write to the query and reference lists
                self._maybe_write_list(ql, path_qry)
                self._maybe_write_list(rl, path_ref)

                # Run FastANI
                result = self.run_proc(q, r, path_qry, path_ref, path_out)
                q_results.put((job, result), timeout=MP_TIMEOUT)
                q_writer.put(True, timeout=MP_TIMEOUT)
            finally:
                shutil.rmtree(dir_tmp)
        return True

    def _writer(self, q_writer, n_total):
        """The writer function, which reports the progress of the workers.

        Parameters
        ----------
        q_writer : multiprocessing.Queue
            A queue of genome ids which have been processed.
        n_total : int
            The total number of items to be processed.
        """
        processed_items = 0
        while True:
            result = q_writer.get(block=True, timeout=MP_TIMEOUT)
            if result is None:
                break
            processed_items += 1
            status = f'==> Processing {processed_items} of {n_total} ' \
                     f'({float(processed_items) * 100 / n_total:.1f}%) comparisons.'
            sys.stdout.write('\r%s' % status)
            sys.stdout.flush()
        sys.stdout.write('\n')
        return True

    def run_proc(self, q, r, ql, rl, output):
        """Runs the FastANI process.

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
        Dict[str, Dict[str, float]]
            The ANI/AF of the query genomes to the reference genomes.
        """
        args = ['fastANI']
        if q is not None:
            args.extend(['-q', q])
        if r is not None:
            args.extend(['-r', r])
        if ql is not None:
            args.extend(['--ql', ql])
        if rl is not None:
            args.extend(['--rl', rl])
        args.extend(['-o', output])

        proc = subprocess.Popen(args, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE, encoding='utf-8')
        stdout, stderr = proc.communicate()

        if proc.returncode != 0:
            self.logger.error('STDOUT:\n' + stdout)
            self.logger.error('STDERR:\n' + stderr)
            raise GTDBTkExit('FastANI returned a non-zero exit code.')

        # Parse the output
        return self._parse_output_file(output)

    def _parse_output_file(self, path_out):
        """Parses the resulting output file from FastANI.

        Parameters
        ----------
        path_out : str
            The path where the output file resides.

        Returns
        -------
        Dict[str, Dict[str, Tuple[float, float]]]
            The ANI/AF of the query genomes to the reference genomes.
        """
        out = dict()
        if os.path.isfile(path_out):
            with open(path_out, 'r') as fh:
                for line in fh.readlines():
                    path_qry, path_ref, ani, frac1, frac2 = line.strip().split('\t')
                    af = round(float(frac1) / float(frac2), 2)
                    if path_qry not in out:
                        out[path_qry] = {path_ref: (float(ani), af)}
                    elif path_ref not in out[path_qry]:
                        out[path_qry][path_ref] = (float(ani), af)
        return out

    def _maybe_write_list(self, d_genomes, path):
        """Writes a query/reference list to disk.

        Parameters
        ----------
        d_genomes : Dict[str, str]
            A dictionary containing the key as the accession and value as path.
        path : str
            The path to write the file to.
        """
        if d_genomes is None or path is None:
            return
        with open(path, 'w') as fh:
            for gid, gid_path in d_genomes.items():
                fh.write(f'{gid_path}\n')

    def _parse_result_queue(self, q_results, path_to_gid):
        """Creates the output dictionary given the results from FastANI

        Parameters
        ----------
        q_results : Queue
            A multiprocessing queue containing raw results.
        path_to_gid : Dict[str ,str]
            A dictionary containing the file path to genome id.

        Returns
        -------
        Dict[str, Dict[str, Dict[str, float]]]
            The ANI/AF of the query genome to all reference genomes.
        """
        out = dict()
        while True:
            q_item = q_results.get(block=True, timeout=MP_TIMEOUT)
            if q_item is None:
                break

            job, result = q_item
            qry_gid = job['qry']

            for path_a, dict_b in result.items():
                for path_b, (ani, af) in dict_b.items():
                    gid_a, gid_b = path_to_gid[path_a], path_to_gid[path_b]

                    # This was done in the forward direction.
                    if gid_a == qry_gid:
                        ref_gid = gid_b
                    # This was done in the reverse direction.
                    elif gid_b == qry_gid:
                        ref_gid = gid_a
                    else:
                        raise GTDBTkExit('FastANI results are malformed.')

                    # Take the largest ANI / AF from either pass.
                    if qry_gid not in out:
                        out[qry_gid] = {ref_gid: {'ani': ani, 'af': af}}
                    elif ref_gid not in out[qry_gid]:
                        out[qry_gid][ref_gid] = {'ani': ani, 'af': af}
                    else:
                        out[qry_gid][ref_gid]['ani'] = max(out[qry_gid][ref_gid]['ani'], ani)
                        out[qry_gid][ref_gid]['af'] = max(out[qry_gid][ref_gid]['af'], af)

        return out
