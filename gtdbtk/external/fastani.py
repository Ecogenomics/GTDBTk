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
        self.version = self._get_version()
        self.minFrac = self._isMinFrac_present()
        self._suppress_v1_warning = False

    @staticmethod
    def _get_version():
        """Returns the version of FastANI on the system path.

        Returns
        -------
        str
            The string containing the fastANI version.
        """
        try:
            proc = subprocess.Popen(['fastANI', '-v'], stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE, encoding='utf-8')
            stdout, stderr = proc.communicate()
            if stderr.startswith('Unknown option:'):
                return 'unknown (<1.3)'
            version = re.search(r'version (.+)', stderr)
            return version.group(1)
        except Exception:
            return 'unknown'
     
    @staticmethod   
    def _isMinFrac_present():
        """Returns true if --minFraction is an option of FastANI on the system path.

        Returns
        -------
        bool
            True/False.
        """
        try:
            proc = subprocess.Popen(['fastANI', '-h'], stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE, encoding='utf-8')
            stdout, stderr = proc.communicate()
            if '--minFraction' in stderr:
                return True
            return False
        except Exception:
            return False

    def run(self, dict_compare, dict_paths):
        """Runs FastANI in batch mode.

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

                    q_worker.put(fwd_dict)
                    q_worker.put(rev_dict)
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

                q_worker.put(fwd_dict)
                q_worker.put(rev_dict)
                n_total += 2

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
                    raise GTDBTkExit('FastANI returned a non-zero exit code.')

            # Stop the writer thread.
            q_writer.put(None)
            p_writer.join()

        except Exception:
            for p in p_workers:
                p.terminate()
            p_writer.terminate()
            raise

        # Process and return each of the results obtained
        path_to_gid = {v: k for k, v in dict_paths.items()}
        q_results.put(None)
        return self._parse_result_queue(q_results, path_to_gid)

    def _worker(self, q_worker, q_writer, q_results):
        """Operates FastANI in list mode.

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
            q = list(job.get('q').values())[0] if job.get('q') is not None else None
            r = list(job.get('r').values())[0] if job.get('r') is not None else None
            ql = job.get('ql')
            rl = job.get('rl')

            # Create a temporary directory to write the lists to.
            with tempfile.TemporaryDirectory(prefix='gtdbtk_fastani_tmp') as dir_tmp:
                path_qry = os.path.join(dir_tmp, 'ql.txt') if ql is not None else None
                path_ref = os.path.join(dir_tmp, 'rl.txt') if rl is not None else None
                path_out = os.path.join(dir_tmp, 'output.txt')

                # Write to the query and reference lists
                self._maybe_write_list(ql, path_qry)
                self._maybe_write_list(rl, path_ref)

                # Run FastANI
                result = self.run_proc(q, r, path_qry, path_ref, path_out)
                q_results.put((job, result))
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
        dict[str, dict[str, float]]
            The ANI/AF of the query genomes to the reference genomes.
        """
        args = ['fastANI']
        if self.minFrac:
            args.extend(['--minFraction', '0'])
        if q is not None:
            args.extend(['-q', q])
        if r is not None:
            args.extend(['-r', r])
        if ql is not None:
            args.extend(['--ql', ql])
        if rl is not None:
            args.extend(['--rl', rl])
        args.extend(['-o', output])
        #self.logger.debug(' '.join(args))
        proc = subprocess.Popen(args, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE, encoding='utf-8')
        stdout, stderr = proc.communicate()

        if proc.returncode != 0:
            self.logger.error('STDOUT:\n' + stdout)
            self.logger.error('STDERR:\n' + stderr)
            raise GTDBTkExit('FastANI returned a non-zero exit code.')

        # Parse the output
        return self.parse_output_file(output)

    def parse_output_file(self, path_out):
        """Parses the resulting output file from FastANI.

        Parameters
        ----------
        path_out : str
            The path where the output file resides.

        Returns
        -------
        dict[str, dict[str, tuple[float, float]]]
            The ANI/AF of the query genomes to the reference genomes.
        """
        out = dict()
        if os.path.isfile(path_out):
            with open(path_out, 'r') as fh:
                for line in fh.readlines():
                    """FastANI version >=1.1 uses tabs instead of spaces to separate columns.
                    Preferentially try split with tabs first instead of split() in-case of 
                    spaces in the file path."""
                    try:
                        try:
                            path_qry, path_ref, ani, frac1, frac2 = line.strip().split('\t')
                        except ValueError:
                            path_qry, path_ref, ani, frac1, frac2 = line.strip().split(' ')
                            if not self._suppress_v1_warning:
                                self.logger.warning('You are using FastANI v1.0, it is recommended '
                                                    'that you update to a more recent version.')
                                self._suppress_v1_warning = True
                        af = float(frac1) / float(frac2)
                        if path_qry not in out:
                            out[path_qry] = {path_ref: (float(ani), af)}
                        elif path_ref not in out[path_qry]:
                            out[path_qry][path_ref] = (float(ani), af)
                    except Exception as e:
                        self.logger.error(f'Exception reading FastANI output: {repr(e)}')
                        raise GTDBTkExit(f'Unable to read line "{line}"')
        return out

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

    def _parse_result_queue(self, q_results, path_to_gid):
        """Creates the output dictionary given the results from FastANI

        Parameters
        ----------
        q_results : mp.Queue
            A multiprocessing queue containing raw results.
        path_to_gid : dict[str ,str]
            A dictionary containing the file path to genome id.

        Returns
        -------
        dict[str, dict[str, dict[str, float]]]
            The ANI/AF of the query genome to all reference genomes.
        """
        out = dict()
        while True:
            q_item = q_results.get(block=True)
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
                        out[qry_gid][ref_gid]['ani'] = max(
                            out[qry_gid][ref_gid]['ani'], ani)
                        out[qry_gid][ref_gid]['af'] = max(
                            out[qry_gid][ref_gid]['af'], af)

        return out
