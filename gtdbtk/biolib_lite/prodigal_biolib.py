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

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2014'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'

import logging
import ntpath
import os
import shutil
import tempfile
from collections import defaultdict, namedtuple

import numpy as np

from .common import remove_extension, make_sure_path_exists, check_file_exists
from .execute import check_on_path
from .parallel import Parallel
from .seq_io import read_fasta, write_fasta


class Prodigal(object):
    """Wrapper for running Prodigal in parallel."""

    def __init__(self, cpus, verbose=True):
        """Initialization.

        Parameters
        ----------
        cpus : int
            Number of cpus to use.
        verbose : boolean
            Flag indicating if progress should be reported.
        """

        self.logger = logging.getLogger('timestamp')

        check_on_path('prodigal')

        self.cpus = cpus
        self.verbose = verbose

    def _producer(self, genome_file):
        """Apply prodigal to genome with most suitable translation table.

        Parameters
        ----------
        genome_file : str
            Fasta file for genome.
        """

        genome_id = remove_extension(genome_file)

        aa_gene_file = os.path.join(self.output_dir, genome_id + '_genes.faa')
        nt_gene_file = os.path.join(self.output_dir, genome_id + '_genes.fna')
        gff_file = os.path.join(self.output_dir, genome_id + '.gff')

        best_translation_table = -1
        table_coding_density = {4: -1, 11: -1}
        if self.called_genes:
            os.system('cp %s %s' %
                      (os.path.abspath(genome_file), aa_gene_file))
        else:
            seqs = read_fasta(genome_file)

            if len(seqs) == 0:
                self.logger.warning('Cannot call Prodigal on an empty genome. '
                                    'Skipped: {}'.format(genome_file))
                return None

            with tempfile.TemporaryDirectory('gtdbtk_prodigal_tmp_') as tmp_dir:

                # determine number of bases
                total_bases = 0
                for seq in seqs.values():
                    total_bases += len(seq)

                # call genes under different translation tables
                translation_tables = [4, 11]

                for translation_table in translation_tables:
                    os.makedirs(os.path.join(tmp_dir, str(translation_table)))
                    aa_gene_file_tmp = os.path.join(tmp_dir, str(
                        translation_table), genome_id + '_genes.faa')
                    nt_gene_file_tmp = os.path.join(tmp_dir, str(
                        translation_table), genome_id + '_genes.fna')
                    gff_file_tmp = os.path.join(tmp_dir, str(
                        translation_table), genome_id + '.gff')

                    # check if there is sufficient bases to calculate prodigal
                    # parameters
                    if total_bases < 100000 or self.meta:
                        proc_str = 'meta'  # use best pre-calculated parameters
                    else:
                        proc_str = 'single'  # estimate parameters from data

                    # If this is a gzipped genome, re-write the uncompressed genome
                    # file to disk
                    prodigal_input = genome_file
                    if genome_file.endswith('.gz'):
                        prodigal_input = os.path.join(
                            tmp_dir, os.path.basename(genome_file[0:-3]) + '.fna')
                        write_fasta(seqs, prodigal_input)

                    # there may be ^M character in the input file,
                    # the following code is similar to dos2unix command to remove
                    # those special characters.
                    with open(prodigal_input, 'r') as fh:
                        text = fh.read().replace('\r\n', '\n')
                    processed_prodigal_input = os.path.join(
                        tmp_dir, os.path.basename(prodigal_input))
                    with open(processed_prodigal_input, 'w') as fh:
                        fh.write(text)

                    args = '-m'
                    if self.closed_ends:
                        args += ' -c'

                    cmd = 'prodigal %s -p %s -q -f gff -g %d -a %s -d %s -i %s > %s 2> /dev/null' % (args,
                                                                                                     proc_str,
                                                                                                     translation_table,
                                                                                                     aa_gene_file_tmp,
                                                                                                     nt_gene_file_tmp,
                                                                                                     processed_prodigal_input,
                                                                                                     gff_file_tmp)
                    os.system(cmd)

                    # determine coding density
                    prodigalParser = ProdigalGeneFeatureParser(gff_file_tmp)

                    codingBases = 0
                    for seq_id, _seq in seqs.items():
                        codingBases += prodigalParser.coding_bases(seq_id)

                    codingDensity = float(codingBases) / total_bases
                    table_coding_density[translation_table] = codingDensity

                # determine best translation table
                if not self.translation_table:
                    best_translation_table = 11
                    if (table_coding_density[4] - table_coding_density[11] > 0.05) and table_coding_density[4] > 0.7:
                        best_translation_table = 4
                else:
                    best_translation_table = self.translation_table

                shutil.copyfile(os.path.join(tmp_dir, str(best_translation_table),
                                             genome_id + '_genes.faa'), aa_gene_file)
                shutil.copyfile(os.path.join(tmp_dir, str(best_translation_table),
                                             genome_id + '_genes.fna'), nt_gene_file)
                shutil.copyfile(os.path.join(tmp_dir, str(best_translation_table),
                                             genome_id + '.gff'), gff_file)

        return (genome_id, aa_gene_file, nt_gene_file, gff_file, best_translation_table, table_coding_density[4],
                table_coding_density[11])

    def _consumer(self, produced_data, consumer_data):
        """Consume results from producer processes.

         Parameters
        ----------
        produced_data : tuple
            Summary statistics for called genes for a specific genome.
        consumer_data : list
            Summary statistics of called genes for each genome.

        Returns
        -------
        consumer_data: d[genome_id] -> namedtuple(aa_gene_file,
                                                    nt_gene_file,
                                                    gff_file,
                                                    best_translation_table,
                                                    coding_density_4,
                                                    coding_density_11)
            Summary statistics of called genes for each genome.
        """

        ConsumerData = namedtuple(
            'ConsumerData',
            'aa_gene_file nt_gene_file gff_file best_translation_table coding_density_4 coding_density_11')
        if consumer_data is None:
            consumer_data = defaultdict(ConsumerData)

        genome_id, aa_gene_file, nt_gene_file, gff_file, best_translation_table, coding_density_4, coding_density_11 = produced_data

        consumer_data[genome_id] = ConsumerData(aa_gene_file,
                                                nt_gene_file,
                                                gff_file,
                                                best_translation_table,
                                                coding_density_4,
                                                coding_density_11)

        return consumer_data

    def _progress(self, processed_items, total_items):
        """Report progress of consumer processes.

        Parameters
        ----------
        processed_items : int
            Number of genomes processed.
        total_items : int
            Total number of genomes to process.

        Returns
        -------
        str
            String indicating progress of data processing.
        """

        return self.progress_str % (processed_items, total_items, float(processed_items) * 100 / total_items)

    def run(self,
            genome_files,
            output_dir,
            called_genes=False,
            translation_table=None,
            meta=False,
            closed_ends=False):
        """Call genes with Prodigal.

        Call genes with prodigal and store the results in the
        specified output directory. For convenience, the
        called_gene flag can be used to indicate genes have
        previously been called and simply need to be copied
        to the specified output directory.

        Parameters
        ----------
        genome_files : list of str
            Nucleotide fasta files to call genes on.
        called_genes : boolean
            Flag indicating if genes are already called.
        translation_table : int
            Specifies desired translation table, use None to automatically
            select between tables 4 and 11.
        meta : boolean
            Flag indicating if prodigal should call genes with the metagenomics procedure.
        closed_ends : boolean
            If True, do not allow genes to run off edges (throws -c flag).
        output_dir : str
            Directory to store called genes.

        Returns
        -------
        d[genome_id] -> namedtuple(best_translation_table
                                            coding_density_4
                                            coding_density_11)
            Summary statistics of called genes for each genome.
        """

        self.called_genes = called_genes
        self.translation_table = translation_table
        self.meta = meta
        self.closed_ends = closed_ends
        self.output_dir = output_dir

        make_sure_path_exists(self.output_dir)

        progress_func = None
        if self.verbose:
            file_type = 'genomes'
            self.progress_str = '  Finished processing %d of %d (%.2f%%) genomes.'
            if meta:
                file_type = 'scaffolds'
                if len(genome_files):
                    file_type = ntpath.basename(genome_files[0])

                self.progress_str = '  Finished processing %d of %d (%.2f%%) files.'

            self.logger.info('Identifying genes within %s: ' % file_type)
            progress_func = self._progress

        parallel = Parallel(self.cpus)
        summary_stats = parallel.run(
            self._producer, self._consumer, genome_files, progress_func)

        # An error was encountered during Prodigal processing, clean up.
        if not summary_stats:
            shutil.rmtree(self.output_dir)

        return summary_stats


class ProdigalGeneFeatureParser(object):
    """Parses prodigal gene feature files (GFF) output."""

    def __init__(self, filename):
        """Initialization.

        Parameters
        ----------
        filename : str
            GFF file to parse.
        """
        check_file_exists(filename)

        self.genes = {}
        self.last_coding_base = {}

        self.__parseGFF(filename)

        self.coding_base_masks = {}
        for seq_id in self.genes:
            self.coding_base_masks[seq_id] = self.__build_coding_base_mask(
                seq_id)

    def __parseGFF(self, filename):
        """Parse genes from GFF file.

        Parameters
        ----------
        filename : str
            GFF file to parse.
        """
        bGetTranslationTable = True
        with open(filename, 'r') as fh:
            for line in fh:
                if bGetTranslationTable and line.startswith('# Model Data'):
                    data_model_info = line.split(':')[1].strip().split(';')
                    dict_data_model = {}
                    for item in data_model_info:
                        k = item.split('=')[0]
                        v = item.split('=')[1]
                        dict_data_model[k] = v

                    self.translationTable = int(
                        dict_data_model.get('transl_table'))
                    bGetTranslationTable = False

                if line[0] == '#':
                    continue

                line_split = line.split('\t')
                seq_id = line_split[0]
                if seq_id not in self.genes:
                    geneCounter = 0
                    self.genes[seq_id] = {}
                    self.last_coding_base[seq_id] = 0

                geneId = seq_id + '_' + str(geneCounter)
                geneCounter += 1
                start = int(line_split[3])
                end = int(line_split[4])

                self.genes[seq_id][geneId] = [start, end]
                self.last_coding_base[seq_id] = max(
                    self.last_coding_base[seq_id], end)

    def __build_coding_base_mask(self, seq_id):
        """Build mask indicating which bases in a sequences are coding.

        Parameters
        ----------
        seq_id : str
            Unique id of sequence.
        """

        # safe way to calculate coding bases as it accounts
        # for the potential of overlapping genes
        coding_base_mask = np.zeros(self.last_coding_base[seq_id], dtype=np.bool)
        for pos in self.genes[seq_id].values():
            coding_base_mask[pos[0]:pos[1] + 1] = True

        return coding_base_mask

    def coding_bases(self, seq_id, start=0, end=None):
        """Calculate number of coding bases in sequence between [start, end).

        To process the entire sequence set start to 0, and
        end to None.

        Parameters
        ----------
        seq_id : str
            Unique id of sequence.
        start : int
            Start calculation at this position in sequence.
        end : int
            End calculation just before this position in the sequence.
        """

        # check if sequence has any genes
        if seq_id not in self.genes:
            return 0

        # set end to last coding base if not specified
        if end is None:
            end = self.last_coding_base[seq_id]

        return np.sum(self.coding_base_masks[seq_id][start:end])
