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
import datetime
import logging
import multiprocessing as mp
import os
import re
import shutil
import subprocess
import tempfile
import time
from pathlib import Path

from tqdm import tqdm

from gtdbtk.exceptions import GTDBTkExit
from gtdbtk.tools import tqdm_log
from gtdbtk.config.common import CONFIG

from collections import namedtuple

SkaniResult = namedtuple('SkaniResult', ['ref_file', 'query_file', 'ani', 'af_r', 'af_q'])


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

    def parse_result_line(self, line):
        """Parses a line from Skani output into a SkaniResult."""
        tokens = line.strip().split('\t')
        if len(tokens) < 5:
            raise ValueError(f"Unexpected Skani result format: {line}")
        return SkaniResult(
            ref_file=tokens[0],
            query_file=tokens[1],
            ani=float(tokens[2]),
            af_r=float(tokens[3]),
            af_q=float(tokens[4])
        )


    def run_vs_all_reps(self, genomes, ref_genomes, prefix,skani_preset=None, report_progress=True):
        """Runs skani against all representatives in the GTDB.

        Parameters
        ----------
        genomes : dict[str, str]
            A dictionary containing the genome ids and their paths.
        ref_genomes : dict[str, str]
            A dictionary containing the reference genome ids and their paths.
        prefix : str
            The prefix to use for the output files.
        preset : str, optional
            The preset to use for skani, e.g. '--medium' or '--slow'.
        report_progress : bool, optional
            If True, report the progress of the skani run.

        Returns
        -------
        dict[str, dict[str, dict[str, float]]]
            A dictionary containing the ANI and AF for each comparison.
        """

        # debug only tke 5 genomes and 5 ref_genomes
        genomes = {k: genomes[k] for k in list(genomes)}
        ref_genomes = {k: ref_genomes[k] for k in list(ref_genomes)}

        # genomes = ref_genomes

        # Create a temporary directory to store the query and reference lists.
        with tempfile.TemporaryDirectory(prefix=prefix) as tmpdir:
            # Create the query and reference lists.
            ql = os.path.join(tmpdir, 'query_list.txt')
            rl = os.path.join(tmpdir, 'ref_list.txt')

            # Write the query and reference lists to disk.
            reverse_dict_ql=self.write_list_to_file(genomes, ql)
            reverse_dict_rl=self.write_list_to_file(ref_genomes, rl)

            # Run skani
            results_all_vs_all= self.run_all_vs_all(ql, rl,reverse_dict_ql,reverse_dict_rl,skani_preset,
                            report_progress=report_progress)

            return self._parse_results(iter(results_all_vs_all))






    def run_all_vs_all(self, ql, rl,reverse_ql,reverse_rl,skani_preset=None, report_progress=True):
        """Runs skani in batch mode against all genomes in the GTDB.

        Parameters
        ----------
        ql : str
            The path to the query list file.
        rl : str
            The path to the reference list file.
        preset : str, optional
            The preset to use for skani, e.g. '--medium' or '--slow'.
        report_progress : bool, optional
            If True, report the progress of the skani run.

        Returns
        -------
        set[tuple[str, str, float, float, float]]
            A set containing tuples of the form (query_id, reference_id, ANI, AF_reference, AF_query).

        """

        ani_af = set()

        #
        args = ['skani', 'dist']
        if skani_preset:
            # is skani_preset doesnt start with "--" then add "--" to it
            if not skani_preset.startswith('--'):
                preset = f'--{skani_preset}'
                args.append(preset)


        #args += ['-t',f'{self.cpus}','-s',f'{skani_s}','--min-af',f'{skani_min_af}','--trace','--ql', ql, '--rl', rl, '-o', '/dev/stdout']
        args += ['-t',f'{self.cpus}','--trace','--ql', ql, '--rl', rl, '-o', '/dev/stdout']

        # print the type of args for debugging
        self.logger.debug('Running skani with the following arguments:')
        self.logger.debug(f'args: {type(args)}')
        # print the args for debugging


        # Read file lists into sets for fast lookup
        with open(ql) as f:
            query_files = set(Path(line.strip()).name for line in f if line.strip())
        with open(rl) as f:
            ref_files = set(Path(line.strip()).name for line in f if line.strip())

        query_total = len(query_files)
        ref_total = len(ref_files)

        # Setup tqdm progress bars
        self.logger.info(f'Sketching genomes')
        ref_bar = tqdm(
            total=ref_total,
            desc="Sketching references",
            unit="genome",
            ncols=100,
            leave=False
        )
        query_bar = tqdm(
            total=query_total,
            desc="Sketching queries",
            unit="genome",
            ncols=100,
            leave=False
        )
        count_ratio = 0



        # Track what's already been counted (avoid duplicate progress)
        seen_queries = set()
        seen_refs = set()
        result_lines = set()
        capture_line=False
        sketch_done = False
        capture_output= []

        proc = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, encoding='utf-8')

        comparison_bar = None

        for line in proc.stdout:
            line = line.strip()
            capture_output.append(line)
            self.logger.debug(line)

            if "TRACE" in line and "Sketching" in line:
                # Extract the file name from path
                try:
                    file_path = line.split()[-2].strip()
                    filename = Path(file_path).name

                    if filename in query_files and filename not in seen_queries:
                        seen_queries.add(filename)
                        query_bar.update(1)

                    elif filename in ref_files and filename not in seen_refs:
                        seen_refs.add(filename)
                        ref_bar.update(1)
                except Exception as e:
                    self.logger.warning(f"Could not parse line: {line} ({e})")

            elif "Generating sketch time" in line:
                time_info = line.split()[-1].strip()
                #convert time_info to seconds if needed
                time_in_seconds = float(time_info)
                # Convert to HH:MM:SS
                minutes = int(time_in_seconds // 60)
                seconds = int(time_in_seconds % 60)

                self.logger.info(f"Sketches done: {minutes}min {seconds}secs")

                query_bar.clear()
                query_bar.close()
                ref_bar.clear()
                ref_bar.close()
                self.logger.info(f"Running comparisons")
                sketch_done = True

                # after sketching is done
                comparison_bar = tqdm(
                    desc="Running comparisons",
                    unit="comparisons",
                    ncols=100,
                    leave=False,  # keep final summary line
                    total=None  # unknown total, just count
                )


            elif sketch_done and ((" TRACE " in line and " Ratio " in line) or
                                  (" DEBUG " in line and " Ref_file " in line and " Query_file " in line)):
                if comparison_bar is not None:
                    comparison_bar.update(1)


            # all Ref_file,Query_file,ANI,Align_fraction_ref,Align_fraction_query,Ref_name,Query_name in line
            elif all(token in line for token in ['Ref_file', 'Query_file', 'ANI', 'Align_fraction_ref', 'Align_fraction_query']):
                # This is the header line, skip it but from now on, we capture the output
                capture_line = True
                if comparison_bar is not None:
                    comparison_bar.close()
                self.logger.info(f'Comparisons finished, capturing results.')
                continue
            elif capture_line and line and "ANI calculation time:" not in line:
                result_lines.add((line.strip()))

        proc.wait()
        #remove the last printed line
  # overwrite with spaces and return to start

        if proc.returncode != 0:
            self.logger.error('STDOUT:\n' + "\n".join(capture_output))
            raise GTDBTkExit('skani returned a non-zero exit code.')

        if len(result_lines)> 1:
            for result_line in result_lines:
                try:
                    parsed = self.parse_result_line(result_line)
                    ani_af.add((
                        reverse_ql.get(parsed.query_file),
                        reverse_rl.get(parsed.ref_file),
                        parsed.ani,
                        parsed.af_r,
                        parsed.af_q
                    ))
                except Exception as e:
                    self.logger.error(f'Error parsing line: {result_line} ({e})')
                    continue
        return ani_af


    def run(self, dict_compare, dict_paths,skani_preset=None, report_progress=True):
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
            q_worker.put((qry, ref, dict_paths.get(qry), dict_paths.get(ref), skani_preset))

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
        return self._parse_results(q_results, from_queue=True)

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
            result = self.run_proc(q, r, q_path, r_path,preset)
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

    def run_proc(self, qid, rid, ql, rl,skani_preset, report_progress=True):
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
        #
        args = ['skani', 'dist']
        if skani_preset:
            # is skani_preset doesnt start with "--" then add "--" to it
            if not skani_preset.startswith('--'):
                preset = f'--{skani_preset}'
                args.append(preset)
        #args += ['-s',f'{skani_s}','--min-af',f'{skani_min_af}','-q', ql, '-r', rl, '-o', '/dev/stdout']
        args += ['-q', ql, '-r', rl, '-o', '/dev/stdout']


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
            parsed = self.parse_result_line(result_lines[1])
            ani_af = (qid, rid, parsed.ani, parsed.af_r, parsed.af_q)
        else:
        #     # genomes too divergent to determine ANI and AF
        #     # with skani so default to zeros
            ani_af = "null"
        #     ani_af = (qid, rid, 0.0, 0.0, 0.0)


        return ani_af


    def write_list_to_file(self, d_genomes, path):
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
        reverse_dict_ql = {v: k for k, v in d_genomes.items()}
        with open(path, 'w') as fh:
            for gid, gid_path in d_genomes.items():
                fh.write(f'{gid_path}\n')
        return reverse_dict_ql

    def _parse_results(self, results, from_queue=False):
        """Parses results into a nested dictionary of ANI and AF values.

        Parameters
        ----------
        results : Iterable or multiprocessing.Queue
            The input results: either an iterable of tuples or a multiprocessing.Queue.
        from_queue : bool
            If True, the input is treated as a Queue and will read until a None is received.

        Returns
        -------
        dict[str, dict[str, dict[str, float]]]
            The ANI/AF of the query genome to all reference genomes.
        """
        out = dict()

        while True:
            if from_queue:
                job = results.get(block=True)
                if job == 'null':
                    continue
                if job is None:
                    break
            else:
                try:
                    job = next(results)
                except StopIteration:
                    break

            qid, rid, ani, af_r, af_q = job
            max_af = max(af_r, af_q) / 100

            if qid not in out:
                out[qid] = {rid: {'ani': ani, 'af': max_af}}
            else:
                out[qid][rid] = {'ani': ani, 'af': max_af}

        return out

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

    def _parse_result_set(self, result_set):
        """Parses the result set from skani.

        Parameters
        ----------
        result_set : set[tuple[str, str, float, float, float]]
            The set of results from skani.

        Returns
        -------
        dict[str, dict[str, dict[str, float]]]
            The ANI/AF of the query genome to all reference genomes.
        """
        out = dict()
        for qid, rid, ani, af_r, af_q in result_set:
            max_af = max(af_r, af_q)
            # af is a percent, we need to divide it by 100
            max_af = max_af / 100
            if qid not in out:
                out[qid] = {rid: {'ani': ani, 'af': max_af}}
            else:
                out[qid][rid] = {'ani': ani, 'af': max_af}

        return out
