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
import subprocess
import tempfile

from gtdbtk.biolib_lite.execute import check_dependencies
from gtdbtk.exceptions import GTDBTkException
from gtdbtk.io.marker.copy_number import CopyNumberFileAR122, CopyNumberFileBAC120
from gtdbtk.io.marker.tophit import TopHitPfamFile, TopHitTigrFile
from gtdbtk.tools import tqdm_log


class HmmAligner(object):
    """Runs HMMalign over a set of genomes."""

    def __init__(self,
                 threads,
                 pfam_top_hit_suffix,
                 tigrfam_top_hit_suffix,
                 protein_file_suffix,
                 pfam_hmm_dir,
                 tigrfam_hmm_dir,
                 bac120_markers,
                 ar122_markers):
        """Initialization."""

        check_dependencies(['hmmalign'])

        self.logger = logging.getLogger('timestamp')

        self.threads = threads
        self.pfam_top_hit_suffix = pfam_top_hit_suffix
        self.tigrfam_top_hit_suffix = tigrfam_top_hit_suffix
        self.protein_file_suffix = protein_file_suffix
        self.pfam_hmm_dir = pfam_hmm_dir
        self.tigrfam_hmm_dir = tigrfam_hmm_dir

        self.bac120_markers = bac120_markers
        self.ar122_markers = ar122_markers

        self.marker_path_prefix = {"PFAM": os.path.join(self.pfam_hmm_dir,
                                                        'individual_hmms'),
                                   "TIGRFAM": os.path.join(os.path.dirname(self.tigrfam_hmm_dir),
                                                           'individual_hmms')}

        self.ar122_marker_sizes = None
        self.bac120_marker_sizes = None

        self.version = self.get_version()

    @staticmethod
    def get_version():
        """ get HMMER version."""
        try:
            env = os.environ.copy()
            proc = subprocess.Popen(['hmmalign', '-h'], stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE, env=env, encoding='utf-8')

            output, error = proc.communicate()
            for line in output.split('\n'):
                if line.startswith('# HMMER'):
                    version = line.split(';')[0].replace('# HMMER', '').strip()
                    return version
            return "(version unavailable)"
        except:
            return "(version unavailable)"

    def _get_hmm_sizes(self):
        ar122, bac120 = dict(), dict()
        for marker_d, out_d in ((self.ar122_markers, ar122),
                                (self.bac120_markers, bac120)):
            for marker_type in ('PFAM', 'TIGRFAM'):
                for marker_name in marker_d[marker_type]:
                    marker_path = os.path.join(self.marker_path_prefix[marker_type], marker_name)
                    marker_name_strip = marker_name.replace('.HMM', '').replace('.hmm', '')
                    out_d[marker_name_strip] = self._get_hmm_size(marker_path)
        self.ar122_marker_sizes = ar122
        self.bac120_marker_sizes = bac120
        return

    def align_marker_set(self, db_genome_ids, marker_set_id):
        """Threaded alignment using hmmalign for a given set of genomes.

        Parameters
        ----------
        db_genome_ids : dict
            A dictionary containing the genome ids and aa paths to process.
        marker_set_id : str
            The marker set of these genomes (bac120/ar122).

        Returns
        -------
        dict
            A dictionary of genome_ids -> aligned sequence.
        """

        # Initialise some values for faster processing.
        self._get_hmm_sizes()

        q_worker = mp.Queue()
        q_writer = mp.Queue()

        for gid, gid_dict in db_genome_ids.items():
            q_worker.put((gid, gid_dict.get('aa_gene_path'), marker_set_id))
        [q_worker.put(None) for _ in range(self.threads)]

        manager = mp.Manager()
        out_dict = manager.dict()

        p_workers = [mp.Process(target=self._worker,
                                args=(q_worker, q_writer, out_dict))
                     for _ in range(self.threads)]

        p_writer = mp.Process(target=self._writer,
                              args=(q_writer, len(db_genome_ids)))

        try:
            p_writer.start()
            for p_worker in p_workers:
                p_worker.start()

            for p_worker in p_workers:
                p_worker.join()

                # Gracefully terminate the program.
                if p_worker.exitcode != 0:
                    raise GTDBTkException(
                        'hmmalign returned a non-zero exit code.')

            q_writer.put(None)
            p_writer.join()

        except Exception:
            for p in p_workers:
                p.terminate()

            p_writer.terminate()
            raise

        return {k: v for k, v in out_dict.items()}

    def _worker(self, q_worker, q_writer, out_dict):
        """The worker function, invoked in a process.

        Parameters
        ----------
        q_worker : multiprocessing.Queue
            A queue of tuples containing the genome id, aa path, and marker.
        q_writer : multiprocessing.Queue
            A queue consumed by the writer, to track progress.
        out_dict : dict
            A thread-safe managed dictionary of results.
        """
        job = q_worker.get(block=True, timeout=None)
        while job is not None:
            db_genome_id, aa_gene_path, marker_set_id = job
            out_dict[db_genome_id] = self._run_multi_align(db_genome_id,
                                                           aa_gene_path,
                                                           marker_set_id)
            q_writer.put(db_genome_id)
            job = q_worker.get(block=True, timeout=None)
        return True

    def _writer(self, q_writer, n_genomes):
        """The writer function, which reports the progress of the workers.

        Parameters
        ----------
        q_writer : multiprocessing.Queue
            A queue of genome ids which have been processed.
        n_genomes : int
            The total number of genomes to be processed.
        """
        with tqdm_log(total=n_genomes, unit='genome') as p_bar:
            for _ in iter(q_writer.get, None):
                p_bar.update()

    def _run_multi_align(self, db_genome_id, path, marker_set_id):
        """
        Returns the concatenated marker sequence for a specific genome
        :param db_genome_id: Selected genome
        :param path: Path to the genomic fasta file for the genome
        :param marker_set_id: Unique ID of marker set to use for alignment
        """

        cur_marker_dir = os.path.dirname(os.path.dirname(path))
        pfam_tophit_file = TopHitPfamFile(cur_marker_dir, db_genome_id)
        tigr_tophit_file = TopHitTigrFile(cur_marker_dir, db_genome_id)
        pfam_tophit_file.read()
        tigr_tophit_file.read()

        if marker_set_id == 'bac120':
            copy_number_file = CopyNumberFileBAC120('/dev/null', None)
            marker_size_d = self.bac120_marker_sizes
        elif marker_set_id == 'ar122':
            copy_number_file = CopyNumberFileAR122('/dev/null', None)
            marker_size_d = self.ar122_marker_sizes
        else:
            raise GTDBTkException('Unknown marker set.')

        copy_number_file.add_genome(db_genome_id, path, pfam_tophit_file, tigr_tophit_file)
        single_copy_hits = copy_number_file.get_single_copy_hits(db_genome_id)

        # gather information for all marker genes
        marker_paths = {"PFAM": os.path.join(self.pfam_hmm_dir, 'individual_hmms'),
                        "TIGRFAM": os.path.join(os.path.dirname(self.tigrfam_hmm_dir), 'individual_hmms')}

        marker_dict_original = {}
        if marker_set_id == "bac120":
            for db_marker in sorted(self.bac120_markers):
                marker_dict_original.update({marker.replace(".HMM", "").replace(".hmm", ""):
                                             os.path.join(
                    marker_paths[db_marker], marker)
                    for marker in self.bac120_markers[db_marker]})
        elif marker_set_id == "ar122":
            for db_marker in sorted(self.ar122_markers):
                marker_dict_original.update({marker.replace(".HMM", "").replace(".hmm", ""):
                                             os.path.join(
                    marker_paths[db_marker], marker)
                    for marker in self.ar122_markers[db_marker]})
        elif marker_set_id == "rps23":
            for db_marker in sorted(self.rps23_markers):
                marker_dict_original.update({marker.replace(".HMM", "").replace(".hmm", ""):
                                             os.path.join(
                    marker_paths[db_marker], marker)
                    for marker in self.rps23_markers[db_marker]})

        # Iterate over each of the expected markers and store the gene sequence.
        gene_dict = dict()
        result_align = dict()
        for marker_id, marker_path in marker_dict_original.items():
            hit = single_copy_hits.get(marker_id)
            if hit:
                gene_dict[marker_id] = {"marker_path": marker_path,
                                        "gene": hit['hit'].gene_id,
                                        "gene_seq": hit['seq'],
                                        "bitscore": hit['hit'].bit_score}
            else:
                result_align[marker_id] = '-' * marker_size_d[marker_id]

        # Align the markers.
        result_align.update(self._run_align(gene_dict))

        # we concatenate the aligned markers together and associate them with
        # the genome.
        return ''.join([x[1] for x in sorted(result_align.items())])

    def _run_align(self, marker_dict):
        """
        Run hmmalign for a set of genes for a specific genome. This is run in a temp folder.
        :param marker_dict: list of markers that need to be aligned
        :param genome: specific genome id
        Returns
        --------------
        List of tuple to be inserted in aligned_markers table
        """
        result_genomes_dict = {}
        with tempfile.TemporaryDirectory(prefix='gtdbtk_tmp_') as dir_tmp:
            input_count = 0
            for markerid, marker_info in marker_dict.items():
                hmmalign_gene_input = os.path.join(dir_tmp,
                                                   f'input_gene{input_count}.fa')
                input_count += 1
                with open(hmmalign_gene_input, 'w') as out_fh:
                    out_fh.write(">{0}\n".format(marker_info.get("gene")))
                    out_fh.write("{0}".format(marker_info.get("gene_seq")))
                proc = subprocess.Popen(["hmmalign", "--outformat", "Pfam", marker_info.get(
                    "marker_path"), hmmalign_gene_input], stdout=subprocess.PIPE,
                                        stderr=subprocess.PIPE, encoding='utf-8')
                stdout, stderr = proc.communicate()

                for line in stderr.splitlines():
                    print(line)

                result = self._get_aligned_marker(
                    marker_info.get("gene"), stdout)
                if len(result) < 1:
                    return "TODO"

                result_genomes_dict[markerid] = result

                input_count += 1
        return result_genomes_dict

    def _get_hmm_size(self, path):
        size = 0
        with open(path) as fp:
            for line in fp:
                if line.startswith("LENG  "):
                    size = line.split("  ")[1]
                    break
        return int(size)

    def _get_aligned_marker(self, hit_name, result_file):
        """
        Parse the output of Hmmalign
        :param hit_name: gene name
        :param result_file: output file from Hmmalign
        """
        hit_seq = None
        mask_seq = None

        for line in result_file.splitlines():
            splitline = line.split(" ", 1)
            if splitline[0] == hit_name.split(" ", 1)[0]:
                rsplitline = line.rsplit(" ", 1)
                hit_seq = rsplitline[-1]
            elif line[0:len("#=GC RF")] == "#=GC RF":
                rsplitline = line.rsplit(" ", 1)
                mask_seq = rsplitline[-1]

        if mask_seq is None:
            raise Exception("Unable to get mask from hmm align result file")

        if hit_seq is None:
            return None

        aligned_marker = ''.join([h for h, m in zip(hit_seq, mask_seq) if m == 'x'])
        return aligned_marker
