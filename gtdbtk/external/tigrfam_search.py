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
import sys

from gtdbtk.exceptions import GTDBTkExit
from gtdbtk.biolib_lite.common import make_sure_path_exists
from gtdbtk.tools import sha256, file_has_checksum


class TigrfamSearch(object):
    """Runs TIGRFAM HMMs over a set of genomes."""

    def __init__(self,
                 threads,
                 tigrfam_hmms,
                 protein_file_suffix,
                 tigrfam_suffix,
                 tigrfam_top_hit_suffix,
                 checksum_suffix,
                 output_dir):
        """Initialization."""

        self.cpus_per_genome = 1
        self.threads = threads
        self.tigrfam_hmms = tigrfam_hmms
        self.protein_file_suffix = protein_file_suffix
        self.tigrfam_suffix = tigrfam_suffix
        self.tigrfam_top_hit_suffix = tigrfam_top_hit_suffix
        self.checksum_suffix = checksum_suffix
        self.output_dir = output_dir
        self.logger = logging.getLogger('timestamp')

    def _topHit(self, tigrfam_file):
        """Determine top hits to TIGRFAMs.

        A gene is assigned to a single TIGRFAM
        family. This will be the top hit among
        all TIGRFAM HMMs and pass the threshold
        for the HMM.

        Parameters
        ----------
        tigrfam_file : str
            Name of file containing hits to TIGRFAM HMMs.
        """
        assembly_dir, filename = os.path.split(tigrfam_file)
        genome_id = filename.replace(self.tigrfam_suffix, '')
        output_tophit_file = os.path.join(self.output_dir, genome_id, filename.replace(self.tigrfam_suffix,
                                                                                       self.tigrfam_top_hit_suffix))

        tophits = {}
        with open(tigrfam_file, 'r') as fh_tigrfam:
            for line in fh_tigrfam:
                if line[0] == '#':
                    continue

                line_split = line.split()
                gene_id = line_split[0]
                hmm_id = line_split[3]
                evalue = float(line_split[4])
                bitscore = float(line_split[5])
                if gene_id in tophits:
                    if bitscore > tophits[gene_id][2]:
                        tophits[gene_id] = (hmm_id, evalue, bitscore)
                else:
                    tophits[gene_id] = (hmm_id, evalue, bitscore)

        with open(output_tophit_file, 'w') as fout:
            fout.write('Gene Id\tTop hits (Family id,e-value,bitscore)\n')
            for gene_id, stats in tophits.items():
                hit_str = ','.join(map(str, stats))
                fout.write('%s\t%s\n' % (gene_id, hit_str))

        # calculate checksum
        checksum = sha256(output_tophit_file)
        with open(output_tophit_file + self.checksum_suffix, 'w') as fout:
            fout.write(checksum)

    def _workerThread(self, queueIn, queueOut):
        """Process each data item in parallel."""
        while True:
            queue_next = queueIn.get(block=True, timeout=None)
            if queue_next is None:
                break
            genome_id, gene_file = queue_next

            output_hit_file = os.path.join(self.output_dir, genome_id, '{}{}'.format(genome_id, self.tigrfam_suffix))
            output_tophit_file = os.path.join(self.output_dir, genome_id, '{}{}'.format(genome_id, self.tigrfam_top_hit_suffix))

            # Genome has already been processed
            if file_has_checksum(output_hit_file) and file_has_checksum(output_tophit_file):
                self.logger.info('Skipping result from a previous run: {}'.format(genome_id))

            # Process this genome
            else:
                genome_dir = os.path.join(self.output_dir, genome_id)
                hmmsearch_out = os.path.join(genome_dir, '{}_tigrfam.out'.format(genome_id))
                make_sure_path_exists(genome_dir)

                args = ['hmmsearch', '-o', hmmsearch_out, '--tblout',
                        output_hit_file, '--noali', '--notextw', '--cut_nc',
                        '--cpu', str(self.cpus_per_genome), self.tigrfam_hmms,
                        gene_file]
                proc = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
                proc_out, proc_err = proc.communicate()

                if proc.returncode != 0:
                    queueOut.put((proc.returncode, genome_id, proc_out, proc_err))
                    sys.exit(proc.returncode)

                # calculate checksum
                checksum = sha256(output_hit_file)
                with open(output_hit_file + self.checksum_suffix, 'w') as fout:
                    fout.write(checksum)

                # identify top hit for each gene
                self._topHit(output_hit_file)

            # allow results to be processed or written to file
            queueOut.put((0, genome_id, None, None))

    def _writerThread(self, numDataItems, writerQueue):
        """Store or write results of worker threads in a single thread."""
        processedItems = 0
        while True:
            t = writerQueue.get(block=True, timeout=None)
            if t is None:
                break

            exit_code, genome_id, stdout, stderr = t

            if exit_code != 0:
                proc_msg = '\n' + '=' * 80 + '\n'
                proc_msg += stdout + '\n'
                proc_msg += '-' * 80 + '\n\n' if stderr else ''
                proc_msg += stderr + '\n' if stderr else ''
                proc_msg += '=' * 80 + '\n'
                self.logger.error('hmmsearch returned exit code {} while '
                                  'processing: {}\n{}'.format(exit_code,
                                                              genome_id,
                                                              proc_msg))
            else:
                processedItems += 1
                statusStr = '==> Finished processing %d of %d (%.1f%%) genomes.' % (processedItems,
                                                                                    numDataItems,
                                                                                    float(
                                                                                        processedItems) * 100 / numDataItems)
                sys.stdout.write('\r%s' % statusStr)
                sys.stdout.flush()
        sys.stdout.write('\n')

    def run(self, gene_files):
        """Annotate genes with TIGRFAM HMMs.

        Parameters
        ----------
        gene_files : iterable
            Gene files in FASTA format to process.
        """
        self.cpus_per_genome = max(1, self.threads / len(gene_files))

        # populate worker queue with data to process
        workerQueue = mp.Queue()
        writerQueue = mp.Queue()

        for f in gene_files:
            workerQueue.put(f)

        for _ in range(self.threads):
            workerQueue.put(None)

        try:
            workerProc = [mp.Process(target=self._workerThread, args=(
                workerQueue, writerQueue)) for _ in range(self.threads)]
            writeProc = mp.Process(target=self._writerThread, args=(
                len(gene_files), writerQueue))

            writeProc.start()

            for p in workerProc:
                p.start()

            for p in workerProc:
                p.join()

            writerQueue.put(None)
            writeProc.join()

            for proc in workerProc:
                if proc.exitcode != 0:
                    raise GTDBTkExit('An error was encountered while running hmmsearch.')

        except Exception:
            for p in workerProc:
                p.terminate()

            writeProc.terminate()
            raise
