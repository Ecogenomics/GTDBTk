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
import sys

from gtdbtk.exceptions import ProdigalException
from ..biolib_lite.prodigal_biolib import (Prodigal as BioLibProdigal)


class Prodigal(object):
    """Perform ab initio gene prediction using Prodigal."""

    def __init__(self,
                 threads,
                 proteins,
                 marker_gene_dir,
                 protein_file_suffix,
                 nt_gene_file_suffix,
                 gff_file_suffix,
                 force):
        """Initialize."""

        self.logger = logging.getLogger('timestamp')

        self.threads = threads

        self.proteins = proteins

        self.marker_gene_dir = marker_gene_dir
        self.protein_file_suffix = protein_file_suffix
        self.nt_gene_file_suffix = nt_gene_file_suffix
        self.gff_file_suffix = gff_file_suffix
        self.force = force

    def _run_prodigal(self, genome_id, fasta_path):
        """Run Prodigal.

        Parameters
        ----------
        fasta_path : str
            Path to FASTA file to process.
        :return
            False if an error occurred.
        """

        output_dir = os.path.join(self.marker_gene_dir, genome_id)

        prodigal = BioLibProdigal(1, False)
        summary_stats = prodigal.run(
            [fasta_path], output_dir, called_genes=self.proteins)

        # An error occured in BioLib Prodigal.
        if not summary_stats:
            if self.force:
                return None
            else:
                raise Exception(
                    "An error was encountered while running Prodigal.")

        summary_stats = summary_stats[summary_stats.keys()[0]]

        # rename output files to adhere to GTDB conventions and desired genome
        # ID
        aa_gene_file = os.path.join(
            output_dir, genome_id + self.protein_file_suffix)
        shutil.move(summary_stats.aa_gene_file, aa_gene_file)

        nt_gene_file = None
        gff_file = None
        translation_table_file = None
        if not self.proteins:
            nt_gene_file = os.path.join(
                output_dir, genome_id + self.nt_gene_file_suffix)
            shutil.move(summary_stats.nt_gene_file, nt_gene_file)

            gff_file = os.path.join(
                output_dir, genome_id + self.gff_file_suffix)
            shutil.move(summary_stats.gff_file, gff_file)

            # save translation table information
            translation_table_file = os.path.join(
                output_dir, 'prodigal_translation_table.tsv')
            fout = open(translation_table_file, 'w')
            fout.write('%s\t%d\n' % ('best_translation_table',
                                     summary_stats.best_translation_table))
            fout.write('%s\t%.2f\n' % ('coding_density_4',
                                       summary_stats.coding_density_4 * 100))
            fout.write('%s\t%.2f\n' % ('coding_density_11',
                                       summary_stats.coding_density_11 * 100))
            fout.close()

        return (aa_gene_file, nt_gene_file, gff_file, translation_table_file, summary_stats.best_translation_table)

    def _worker(self, out_dict, worker_queue, writer_queue):
        """This worker function is invoked in a process."""

        while True:
            data = worker_queue.get(block=True, timeout=None)
            if data is None:
                break

            genome_id, file_path = data

            rtn_files = self._run_prodigal(genome_id, file_path)

            # Only proceed if an error didn't occur in BioLib Prodigal
            if rtn_files:
                aa_gene_file, nt_gene_file, gff_file, translation_table_file, best_translation_table = rtn_files
                prodigal_infos = {}
                prodigal_infos["aa_gene_path"] = aa_gene_file
                prodigal_infos["nt_gene_path"] = nt_gene_file
                prodigal_infos["gff_path"] = gff_file
                prodigal_infos["translation_table_path"] = translation_table_file
                prodigal_infos["best_translation_table"] = best_translation_table

                out_dict[genome_id] = prodigal_infos
            writer_queue.put(genome_id)

    def _writer(self, num_items, writer_queue):
        """Store or write results of worker threads in a single thread."""
        processed_items = 0
        while processed_items < num_items:
            a = writer_queue.get(block=True, timeout=None)
            if a is None:
                break

            processed_items += 1
            statusStr = '==> Finished processing %d of %d (%.1f%%) genomes.' % (processed_items,
                                                                                num_items,
                                                                                float(processed_items) * 100 / num_items)
            sys.stdout.write('%s\r' % statusStr)
            sys.stdout.flush()

        sys.stdout.write('\n')

    def run(self, genomic_files):
        """Run Prodigal across a set of genomes.

        Parameters
        ----------
        genomic_files : dict
            Dictionary indicating the genomic and gene file for each genome.
        """

        # populate worker queue with data to process
        worker_queue = mp.Queue()
        writer_queue = mp.Queue()

        for genome_id, file_path in genomic_files.iteritems():
            worker_queue.put([genome_id, file_path])

        for _ in range(self.threads):
            worker_queue.put(None)

        try:
            manager = mp.Manager()
            out_dict = manager.dict()

            worker_proc = [mp.Process(target=self._worker, args=(out_dict,
                                                                 worker_queue,
                                                                 writer_queue)) for _ in range(self.threads)]
            writer_proc = mp.Process(target=self._writer, args=(
                len(genomic_files), writer_queue))

            writer_proc.start()
            for p in worker_proc:
                p.start()

            for p in worker_proc:
                p.join()

                # Gracefully terminate the program.
                if p.exitcode != 0:
                    self.logger.error('Prodigal returned a non-zero exit code.')
                    raise ProdigalException

            writer_queue.put(None)
            writer_proc.join()
        except Exception:
            for p in worker_proc:
                p.terminate()

            writer_proc.terminate()
            self.logger.error('An exception was caught while running Prodigal.')
            raise ProdigalException

        result_dict = {k: v for k, v in out_dict.items()}
        return result_dict
