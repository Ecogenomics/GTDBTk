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

import os
import logging
import tempfile
import shutil
import subprocess
import ntpath
from collections import defaultdict, namedtuple

from gtdbtk.biolib_lite.common import remove_extension, make_sure_path_exists
from gtdbtk.biolib_lite.seq_io import read_fasta, write_fasta
from gtdbtk.biolib_lite.parallel import Parallel
from gtdbtk.biolib_lite.execute import check_on_path
import numpy as np


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
        table_prob = {4: -1, 11: -1}
        if self.called_genes:
            os.system('cp %s %s' % (os.path.abspath(genome_file), aa_gene_file))
        else:
            seqs = read_fasta(genome_file)

            if len(seqs) == 0:
                self.logger.warning('Cannot call Prodigal on an empty genome. Skipped: {}'.format(genome_file))
                return None

            tmp_dir = tempfile.mkdtemp()

            # determine number of bases
            total_bases = 0
            for seq in seqs.values():
                total_bases += len(seq)

            # call genes under different translation tables
            if self.translation_table:
                translation_tables = [self.translation_table]
            else:
                translation_tables = [4, 11]

            translation_table_gffs = dict()
            tln_table_stats = dict()
            for translation_table in translation_tables:
                os.makedirs(os.path.join(tmp_dir, str(translation_table)))
                aa_gene_file_tmp = os.path.join(tmp_dir, str(translation_table), genome_id + '_genes.faa')
                nt_gene_file_tmp = os.path.join(tmp_dir, str(translation_table), genome_id + '_genes.fna')

                # check if there are sufficient bases to calculate prodigal parameters
                if total_bases < 100000 or self.meta:
                    proc_str = 'meta'  # use best precalculated parameters
                else:
                    proc_str = 'single'  # estimate parameters from data

                # If this is a gzipped genome, re-write the uncompressed genome file to disk
                prodigal_input = genome_file
                if genome_file.endswith('.gz'):
                    prodigal_input = os.path.join(tmp_dir, os.path.basename(genome_file[0:-3]) + '.fna')
                    write_fasta(seqs, prodigal_input)

                args = ['prodigal', '-m', '-p', proc_str, '-q', '-f', 'gff', '-g', str(translation_table), '-a',
                        aa_gene_file_tmp, '-d', nt_gene_file_tmp, '-i', prodigal_input]
                if self.closed_ends:
                    args.append('-c')

                proc = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
                proc_out, proc_err = proc.communicate()
                gff_stdout = proc_out

                translation_table_gffs[translation_table] = gff_stdout

                if proc.returncode != 0:
                    self.logger.warning('Prodigal returned a non-zero exit code while processing: {}'.format(genome_file))
                    return None

                # determine coding density
                prodigal_parser = ProdigalGeneFeatureParser(gff_stdout)

                # Skip if no genes were called.
                if prodigal_parser.n_sequences_processed() == 0:
                    shutil.rmtree(tmp_dir)
                    self.logger.warning(
                        'No genes were called! Check the quality of your genome. Skipped: {}'.format(genome_file))
                    return None

                # Save the statistics for this translation table
                prodigal_stats = prodigal_parser.generate_statistics()
                tln_table_stats[translation_table] = prodigal_stats
                table_coding_density[translation_table] = prodigal_stats.coding_density

            # determine best translation table
            if not self.translation_table:

                # Logistic classifier coefficients
                b0 = 12.363017423768538
                bi = np.array(
                    [0.01212327382066545, -0.9250857181041326, -0.10176647009345675, 0.7733711446656522,
                     0.6355731038236031, -0.1631355971443377, -0.14713264317198863, -0.10320909026025472,
                     0.09621494439016824, 0.4992209080695785, 1.159933669041023, -0.0507139271834123,
                     1.2619603455217179, 0.24392226222721214, -0.08567859197118802, -0.18759562346413916,
                     0.13136209122186523, -0.1399459561138417, 2.08086235029142, 0.6917662070950119])

                # Scale x
                scaler_mean = np.array(
                    [0.0027036907781622732, -1.8082140490218692, -8.511942254988097e-08, 19.413811775420918,
                     12.08719100126732, 249.89521467118365, 0.0011868456444391487, -0.0007358432829349235,
                     0.004750880986023392, -0.04096159411654551, -0.12505492579693805, -0.03749033894554058,
                     0.13053986993752234, -0.15914556336256136, -0.6075506034967058, 0.06704648371665446,
                     0.04316693333324335, 0.26905236546875266, 0.010326462563249823, 333.3320678912514])
                scaler_scale = np.array(
                    [0.08442772272873166, 2.043313786484819, 2.917510891467501e-05, 22.577812640992242,
                     12.246767248868036, 368.87834547339907, 0.0014166252200216657, 0.0014582164250905056,
                     0.025127203671053467, 0.5095427815162036, 0.2813128128116135, 0.2559877920464989,
                     1.274371529860827, 0.7314782174742842, 1.6885750374356985, 0.17019369029012987,
                     0.15376309021975043, 0.583965556283342, 0.025076680822882474, 544.3648797867784])
                xi = np.array(tln_table_stats[11]) - np.array(tln_table_stats[4])
                xi -= scaler_mean
                xi /= scaler_scale

                # If xi are all 0, then P(11) = 1.
                prob_tbl_11 = 1 / (1 + np.exp(-1 * (b0 + (bi * xi).sum())))
                best_translation_table = 11 if prob_tbl_11 >= 0.5 else 4
                table_prob[4] = 1.0 - prob_tbl_11
                table_prob[11] = prob_tbl_11

            else:
                best_translation_table = self.translation_table

            shutil.copyfile(os.path.join(tmp_dir, str(
                best_translation_table), genome_id + '_genes.faa'), aa_gene_file)
            shutil.copyfile(os.path.join(tmp_dir, str(
                best_translation_table), genome_id + '_genes.fna'), nt_gene_file)
            with open(gff_file, 'w') as f:
                f.write(translation_table_gffs[best_translation_table])

            # clean up temporary files
            shutil.rmtree(tmp_dir)
        return genome_id, aa_gene_file, nt_gene_file, gff_file, best_translation_table, table_coding_density[4], table_coding_density[11], table_prob[4], table_prob[11]

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
                                                    coding_density_11,
                                                    probability_4,
                                                    probability_11)
            Summary statistics of called genes for each genome.
        """

        ConsumerData = namedtuple('ConsumerData', [
            'aa_gene_file',
            'nt_gene_file',
            'gff_file',
            'best_translation_table',
            'coding_density_4',
            'coding_density_11',
            'probability_4',
            'probability_11'
        ])

        if consumer_data is None:
            consumer_data = defaultdict(ConsumerData)

        genome_id, aa_gene_file, nt_gene_file, gff_file, \
        best_translation_table, coding_density_4, coding_density_11, \
        probability_4, probability_11 = produced_data

        consumer_data[genome_id] = ConsumerData(aa_gene_file,
                                                nt_gene_file,
                                                gff_file,
                                                best_translation_table,
                                                coding_density_4,
                                                coding_density_11,
                                                probability_4,
                                                probability_11)

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

    def __init__(self, file_contents):
        """Initialization.

        Parameters
        ----------
        file_contents : str
            The contents of the GFF file to parse.
        """
        self.sequence_data = dict()
        self.model_data = dict()
        self.called_genes = defaultdict(list)
        self.__parse_gff(file_contents)

        self.coding_masks = dict()
        for seq_id in self.called_genes:
            self.coding_masks[seq_id] = self.__build_coding_mask(seq_id)

    def n_sequences_processed(self):
        """Return how many sequences were processed.

        Parameters
        ----------
        :return : int
            The number of sequences processed.
        """
        return len(self.called_genes)

    def __build_coding_mask(self, seq_id):
        """Construct a numpy array representing the coding mask for a sequence.

        Args:
            seq_id (str): Construct the coding mask for this sequence.

        Returns:
            np.ndarray: True if that position codes for a gene, False otherwise.
        """
        coding_mask = np.zeros(self.sequence_data[seq_id].seqlen, dtype=np.bool_)
        for gene in self.called_genes[seq_id]:
            coding_mask[gene.start - 1:gene.end] = True
        return coding_mask

    def __parse_gff(self, file_contents):
        """Parse genes from GFF file and save them to the member variables.

        Parameters
        ----------
        file_contents : str
            A string containing the GFF contents.
        """
        def parse_sequence_data(cur_line):
            SequenceData = namedtuple('SequenceData', 'seqnum seqlen seqhdr')
            seqnum, seqlen, seqhdr = cur_line.split('# Sequence Data: ')[1].split(';')
            seqnum = int(seqnum.split('=')[1])
            seqlen = int(seqlen.split('=')[1])
            seqhdr = seqhdr.split('=')[1].replace('"', '')
            return SequenceData(seqnum, seqlen, seqhdr)

        def parse_model_data(cur_line):
            ModelData = namedtuple('ModelData', 'version run_type model gc_cont transl_table uses_sd')
            version, run_type, model, gc_cont, transl_table, uses_sd = cur_line.split('# Model Data: ')[1].split(';')
            version = version.split('=')[1]
            run_type = run_type.split('=')[1]
            model = model.split('=')[1].replace('"', '')
            gc_cont = float(gc_cont.split('=')[1])
            transl_table = int(transl_table.split('=')[1])
            uses_sd = int(uses_sd.split('=')[1])
            return ModelData(version, run_type, model, gc_cont, transl_table, uses_sd)

        def parse_gene_data(cur_line):
            GeneData = namedtuple('GeneData', [
                'seqid', 'source', 'type', 'start',  'end',  'score',  'strand', 'phase',
                'id', 'partial', 'start_type', 'rbs_motif', 'rbs_spacer', 'gc_cont', 'conf',
                'cscore', 'sscore', 'rscore', 'uscore', 'tscore'
            ])
            seqid, source, gene_type, start, end, _, strand, phase, attrib = cur_line.split('\t')
            gene_id, partial, start_type, rbs_motif, rbs_spacer, gc_cont, conf, score, cscore, sscore, rscore, uscore, tscore, _ = attrib.split(';')
            start = int(start)
            end = int(end)
            gene_id = gene_id.split('=')[1]
            partial = partial.split('=')[1]
            start_type = start_type.split('=')[1]
            rbs_motif = rbs_motif.split('=')[1]
            rbs_spacer = rbs_spacer.split('=')[1]
            gc_cont = float(gc_cont.split('=')[1])
            conf = float(conf.split('=')[1])
            score = float(score.split('=')[1])
            cscore = float(cscore.split('=')[1])
            sscore = float(sscore.split('=')[1])
            rscore = float(rscore.split('=')[1])
            uscore = float(uscore.split('=')[1])
            tscore = float(tscore.split('=')[1])
            return GeneData(seqid, source, gene_type, start, end, score, strand, phase, gene_id, partial, start_type,
                            rbs_motif, rbs_spacer, gc_cont, conf, cscore, sscore, rscore, uscore, tscore)

        cur_gene_id = None
        for line in file_contents.splitlines():
            if line.startswith('##') or line == '' or line == '"':
                continue

            if line.startswith('# Sequence Data:'):
                cur_sequence_data = parse_sequence_data(line)
                cur_gene_id = cur_sequence_data.seqhdr.split(' ')[0]
                self.sequence_data[cur_gene_id] = cur_sequence_data

            elif line.startswith('# Model Data'):
                cur_model_data = parse_model_data(line)
                self.model_data[cur_gene_id] = cur_model_data

            else:
                cur_gene_data = parse_gene_data(line)
                assert(cur_gene_id == cur_gene_data.seqid)
                self.called_genes[cur_gene_id].append(cur_gene_data)

    def coding_bases(self, seq_id, start=1, end=None):
        """Calculate number of coding bases in sequence between [start, end].

        To process the entire sequence set start to 1, and
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
        if seq_id not in self.called_genes:
            return 0

        # set end to last coding base if not specified
        if end is None:
            end = self.sequence_data[seq_id].seqlen

        return self.coding_masks[seq_id][start-1:end].sum()

    def generate_statistics(self):
        """Generate statistics on the genes called.

        Returns:
            namedtuple: Containing those statistics.
        """
        Stats = namedtuple('Stats', [
            'conf_50', 'conf_std', 'conf_max',
            'cscore_50', 'cscore_std', 'cscore_max',
            'gc_cont_50', 'gc_cont_std', 'gc_cont_max',
            'tscore_50', 'tscore_std', 'tscore_max',
            'rscore_50', 'rscore_std', 'rscore_max',
            'uscore_50', 'uscore_std', 'uscore_max',
            'coding_density', 'genes_called'
        ])

        # Determine coding density
        coding_density_numer = 0
        coding_density_denom = 0
        for seq_id, mask in self.coding_masks.items():
            coding_density_numer += mask.sum()
            coding_density_denom += mask.shape[0]
        coding_density = float(coding_density_numer) / coding_density_denom

        # Determine gene statistics
        genes_called = list()
        gc_cont = list()
        conf = list()
        score = list()
        cscore = list()
        rscore = list()
        uscore = list()
        tscore = list()
        for seq_id, gene_list in self.called_genes.items():
            genes_called.append(len(gene_list))
            for cur_gene in gene_list:
                gc_cont.append(cur_gene.gc_cont)
                conf.append(cur_gene.conf)
                score.append(cur_gene.score)
                cscore.append(cur_gene.cscore)
                rscore.append(cur_gene.rscore)
                uscore.append(cur_gene.uscore)
                tscore.append(cur_gene.tscore)

        get_stat = lambda feature: (np.percentile(feature, 50), np.std(feature), np.max(feature))
        conf_50, conf_std, conf_max = get_stat(conf)
        cscore_50, cscore_std, cscore_max = get_stat(cscore)
        gc_cont_50, gc_cont_std, gc_cont_max = get_stat(gc_cont)
        tscore_50, tscore_std, tscore_max = get_stat(tscore)
        rscore_50, rscore_std, rscore_max = get_stat(rscore)
        uscore_50, uscore_std, uscore_max = get_stat(uscore)
        genes_called = np.sum(genes_called)

        stats = Stats(conf_50, conf_std, conf_max,
                      cscore_50, cscore_std, cscore_max,
                      gc_cont_50, gc_cont_std, gc_cont_max,
                      tscore_50, tscore_std, tscore_max,
                      rscore_50, rscore_std, rscore_max,
                      uscore_50, uscore_std, uscore_max,
                      coding_density, genes_called)

        return stats
