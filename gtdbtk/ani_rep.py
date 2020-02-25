import logging
import os
from collections import defaultdict

from gtdbtk.biolib_lite.execute import check_dependencies
from gtdbtk.config.config import FASTANI_GENOMES, FASTANI_GENOMES_EXT
from gtdbtk.config.output import DIR_ANI_REP_INT_MASH
from gtdbtk.external.fastani import FastANI
from gtdbtk.external.mash import Mash


class ANIRep(object):
    """Computes a list of genomes to a list of representatives."""

    def __init__(self, cpus):
        self.logger = logging.getLogger('timestamp')
        self.cpus = cpus

    @staticmethod
    def check_dependencies(no_mash):
        """Exits the system if the required programs are not on the path."""
        dependencies = ['fastANI']
        if not no_mash:
            dependencies.append('mash')
        check_dependencies(dependencies)

    @staticmethod
    def _get_ref_genomes():
        """Returns a dictionary of genome accession to genome path."""
        ref_genomes = dict()
        for f_name in os.listdir(FASTANI_GENOMES):
            if f_name.endswith(FASTANI_GENOMES_EXT):
                accession = f_name.split(FASTANI_GENOMES_EXT)[0]
                ref_genomes[accession] = os.path.join(FASTANI_GENOMES, f_name)
        return ref_genomes

    def _write_results(self, results, out_dir):
        """Writes the FastANI results to disk."""
        out_path = os.path.join(out_dir, 'fastani_results.tsv')
        with open(out_path, 'w') as fh:
            fh.write('query\treference\tani\taf\n')
            for qry_gid, ref_hits in sorted(results.items()):
                for ref_gid, ref_hit in sorted(ref_hits.items(), key=lambda x: (-x[1]['af'], -x[1]['ani'], x[0])):
                    fh.write(f'{qry_gid}\t{ref_gid}\t{ref_hit["ani"]}\t{ref_hit["af"]}\n')
        self.logger.info(f'Results saved to: {out_path}')

    def run(self, genomes, no_mash, max_d, out_dir, prefix, mash_k, mash_v, mash_s):
        """Runs the pipeline."""
        self.check_dependencies(no_mash)

        ref_genomes = self._get_ref_genomes()
        d_compare = defaultdict(set)
        d_paths = {**genomes, **ref_genomes}

        # Pre-filter using Mash if specified.
        if not no_mash:
            dir_mash = os.path.join(out_dir, DIR_ANI_REP_INT_MASH)
            mash = Mash(self.cpus, dir_mash, prefix)
            self.logger.info(f'Using Mash version {mash.version()}')
            mash_results = mash.run(genomes, ref_genomes, max_d, mash_k, mash_v, mash_s)
            for qry_gid, ref_hits in mash_results.items():
                d_compare[qry_gid] = d_compare[qry_gid].union(set(ref_hits.keys()))

        # Compare against all reference genomes.
        else:
            for qry_gid in genomes:
                d_compare[qry_gid] = set(ref_genomes.keys())

        self.logger.info('Calculating ANI with FastANI.')
        fastani = FastANI(self.cpus, force_single=True)
        fastani_results = fastani.run(d_compare, d_paths)

        self._write_results(fastani_results, out_dir)
