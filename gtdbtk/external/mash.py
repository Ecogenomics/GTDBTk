import logging
import os
import re
import subprocess
import sys
import tempfile
from collections import defaultdict

from gtdbtk.biolib_lite.common import make_sure_path_exists
from gtdbtk.exceptions import GTDBTkExit


class Mash(object):
    """Runs Mash against genomes."""

    def __init__(self, cpus, out_dir, prefix):
        self.logger = logging.getLogger('timestamp')
        self.cpus = max(cpus, 1)
        self.out_dir = out_dir
        self.prefix = prefix

    @staticmethod
    def version():
        """Returns the version of mash, or 'unknown' if not known."""
        try:
            proc = subprocess.Popen(['mash', '--version'], stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE, encoding='utf-8')
            stdout, stderr = proc.communicate()
            if len(stdout) > 0:
                return stdout.strip()
            else:
                return 'unknown'
        except Exception:
            return 'unknown'

    def run(self, qry, ref, mash_d, mash_k, mash_v, mash_s):
        """Run Mash on a set of reference and query genomes.

        Parameters
        ----------
        qry : Dict[str, str]
            The set of query genomes and their path.
        ref : Dict[str, str]
            The set of reference genomes and their path.
        mash_d : float
            The maximum distance to report.
        mash_k : int
            The number of k-mers to store for each sequence.
        mash_v : float
            Maximum p-value to report.
        mash_s: int
            Maximum number of non-redundant hashes.
        """
        qry_sketch = QrySketchFile(qry, self.out_dir, self.prefix, self.cpus, mash_k, mash_s)
        ref_sketch = RefSketchFile(ref, self.out_dir, self.prefix, self.cpus, mash_k, mash_s)

        # Generate an output file comparing the distances between these genomes.
        mash_dists = DistanceFile(qry_sketch, ref_sketch, self.out_dir, self.prefix,
                                  self.cpus, max_d=mash_d, mash_v=mash_v)
        results = mash_dists.read()

        # Convert the results back to the accession
        path_to_qry = {v: k for (k, v) in qry.items()}
        path_to_ref = {v: k for (k, v) in ref.items()}
        out = defaultdict(dict)
        for qry_path, ref_hits in results.items():
            for ref_path, hit in ref_hits.items():
                out[path_to_qry[qry_path]][path_to_ref[ref_path]] = hit
        return out

    def _run(self, ref_msh, qry_msh, max_d):
        args = ['mash', 'dist', '-p', self.cpus, '-d', max_d, ref_msh, qry_msh]
        args = list(map(str, args))
        proc = subprocess.Popen(args, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE, encoding='utf-8')
        stdout, stderr = proc.communicate()
        if proc.returncode != 0:
            raise GTDBTkExit(f'Error running Mash dist: {proc.stderr.read()}')

        out = defaultdict(dict)
        for ref_id, qry_id, dist, p_val, shared_n, shared_d in re.findall(r'(.+)\t(.+)\t(.+)\t(.+)\t(\d+)\/(\d+)\n',
                                                                          stdout):
            dist, p_val = float(dist), float(p_val)
            shared_num, shared_den = int(shared_n), int(shared_d)
            out[qry_id][ref_id] = (dist, p_val, shared_num, shared_den)

        return out


class DistanceFile(object):
    """The resulting distance file from the mash dist command."""
    name = 'mash_distances.tsv'

    def __init__(self, qry_sketch, ref_sketch, root, prefix, cpus, max_d, mash_v):
        self.logger = logging.getLogger('timestamp')
        self.qry_sketch = qry_sketch
        self.ref_sketch = ref_sketch
        self.path = os.path.join(root, f'{prefix}.{self.name}')
        self.cpus = cpus
        self.max_d = max_d
        self.mash_v = mash_v

        self._calculate()

    def _calculate(self):
        self.logger.info('Calculating Mash distances.')
        args = ['mash', 'dist', '-p', self.cpus, '-d', self.max_d, '-v',
                self.mash_v, self.ref_sketch.path, self.qry_sketch.path]
        args = list(map(str, args))
        with open(self.path, 'w') as f_out:
            proc = subprocess.Popen(args, stdout=f_out,
                                    stderr=subprocess.PIPE, encoding='utf-8')
            _, stderr = proc.communicate()
        if proc.returncode != 0:
            raise GTDBTkExit(f'Error running Mash dist: {proc.stderr.read()}')

    def read(self):
        """Reads the results of the distance file.

        Returns
        -------
        Dict[str, Dict[str, Tuple[float, float, int, int]]
        """
        out = defaultdict(dict)
        with open(self.path, 'r') as fh:
            hits = re.findall(r'(.+)\t(.+)\t(.+)\t(.+)\t(\d+)\/(\d+)\n', fh.read())
        for ref_id, qry_id, dist, p_val, shared_n, shared_d in hits:
            dist, p_val = float(dist), float(p_val)
            shared_num, shared_den = int(shared_n), int(shared_d)
            out[qry_id][ref_id] = (dist, p_val, shared_num, shared_den)
        return out


class SketchFile(object):
    """Output files which are generated by mash sketch."""

    def __init__(self, genomes, path, cpus, k, s):
        self.logger = logging.getLogger('timestamp')
        self.genomes = genomes
        self.path = path
        self.data = dict()
        self.args = dict()
        self.cpus = cpus
        self.k = k
        self.s = s

        make_sure_path_exists(os.path.dirname(self.path))

        # Use the pre-existing sketch file, otherwise generate it.
        if os.path.isfile(self.path):
            self.logger.info(f'Loading data from existing Mash sketch file: {self.path}')
            self._load_metadata()
            if not self._is_consistent():
                raise GTDBTkExit(f'The sketch file is not consistent with the '
                                 f'input genomes. Remove the existing sketch '
                                 f'file or specify a new output directory.')
        else:
            self.logger.info(f'Creating Mash sketch file: {self.path}')
            self._generate()

    def _load_metadata(self):
        """Loads the metadata from an existing Mash sketch file."""
        args = ['mash', 'info', '-t', self.path]
        proc = subprocess.Popen(args, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE, encoding='utf-8')
        stdout, stderr = proc.communicate()

        if proc.returncode != 0:
            raise GTDBTkExit(f'Error reading Mash sketch file {self.path}:\n{stderr}')

        for hashes, length, path in re.findall(r'(\d+)\t(\d+)\t(.+)\t.+\n', stdout):
            self.data[path] = (int(hashes), int(length))

    def _is_consistent(self):
        """Returns True if the sketch was generated from the genomes."""
        return set(self.data.keys()) == set(self.genomes.values())

    def _generate(self):
        """Generate a new sketch file."""
        with tempfile.TemporaryDirectory(prefix='gtdbtk_mash_tmp_') as dir_tmp:
            path_genomes = os.path.join(dir_tmp, 'genomes.txt')
            with open(path_genomes, 'w') as fh:
                for path in self.genomes.values():
                    fh.write(f'{path}\n')

            args = ['mash', 'sketch', '-l', '-p', self.cpus, path_genomes, '-o',
                    self.path, '-k', self.k, '-s', self.s]
            args = list(map(str, args))
            proc = subprocess.Popen(args, stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE, encoding='utf-8')
            n_processed = 0
            while True:
                line = proc.stderr.readline()
                if not line:
                    sys.stdout.write('\n')
                    sys.stdout.flush()
                    break
                if line.startswith('Sketching'):
                    n_processed += 1
                    pct = round(100 * (n_processed / len(self.genomes)), 2)
                    sys.stdout.write(f'\r==> Sketching {n_processed} of '
                                     f'{len(self.genomes)} ({pct}%) genomes')
                    sys.stdout.flush()
            proc.wait()

            if proc.returncode != 0 or not os.path.isfile(self.path):
                raise GTDBTkExit(f'Error generating Mash sketch: {proc.stderr.read()}')


class QrySketchFile(SketchFile):
    name = 'user_query_sketch.msh'

    def __init__(self, genomes, root, prefix, cpus, k, s):
        path = os.path.join(root, f'{prefix}.{self.name}')
        super().__init__(genomes, path, cpus, k, s)


class RefSketchFile(SketchFile):
    name = 'gtdb_ref_sketch.msh'

    def __init__(self, genomes, root, prefix, cpus, k, s):
        path = os.path.join(root, f'{prefix}.{self.name}')
        super().__init__(genomes, path, cpus, k, s)
