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
import os
from collections import defaultdict
from shutil import copy
from typing import Dict, Tuple, Optional

import numpy as np

import gtdbtk.config.config as Config
from gtdbtk.biolib_lite.common import make_sure_path_exists
from gtdbtk.biolib_lite.execute import check_dependencies
from gtdbtk.biolib_lite.seq_io import read_fasta
from gtdbtk.biolib_lite.taxonomy import Taxonomy
from gtdbtk.config.output import *
from gtdbtk.exceptions import GenomeMarkerSetUnknown, MSAMaskLengthMismatch, InconsistentGenomeBatch
from gtdbtk.external.pfam_search import PfamSearch
from gtdbtk.external.prodigal import Prodigal
from gtdbtk.external.tigrfam_search import TigrfamSearch
from gtdbtk.io.marker.copy_number import CopyNumberFileAR122, CopyNumberFileBAC120
from gtdbtk.io.marker.tophit import TopHitPfamFile, TopHitTigrFile
from gtdbtk.io.marker_info import MarkerInfoFileAR122, MarkerInfoFileBAC120
from gtdbtk.io.prodigal.tln_table import TlnTableFile
from gtdbtk.io.prodigal.tln_table_summary import TlnTableSummaryFile
from gtdbtk.pipeline import align
from gtdbtk.tools import merge_two_dicts, symlink_f, tqdm_log
from gtdbtk.trim_msa import TrimMSA


class Markers(object):
    """Identify and align marker genes."""

    def __init__(self, cpus=1, debug=False):
        """Initialize."""

        self.logger = logging.getLogger('timestamp')
        self.warnings = logging.getLogger('warnings')

        self.cpus = cpus
        self.debug = debug
        self.marker_gene_dir = None

        self.genome_file_suffix = GENOME_FILE_SUFFIX
        self.protein_file_suffix = PROTEIN_FILE_SUFFIX
        self.nt_gene_file_suffix = NT_GENE_FILE_SUFFIX
        self.gff_file_suffix = GFF_FILE_SUFFIX
        self.checksum_suffix = CHECKSUM_SUFFIX

        self.taxonomy_file = Config.TAXONOMY_FILE

        self.pfam_hmm_dir = Config.PFAM_HMM_DIR
        self.pfam_suffix = PFAM_SUFFIX
        self.pfam_top_hit_suffix = PFAM_TOP_HIT_SUFFIX

        self.tigrfam_hmms = Config.TIGRFAM_HMMS
        self.tigrfam_suffix = TIGRFAM_SUFFIX
        self.tigrfam_top_hit_suffix = TIGRFAM_TOP_HIT_SUFFIX

    def _report_identified_marker_genes(self, gene_dict, outdir, prefix,
                                        write_single_copy_genes):
        """Report statistics for identified marker genes."""

        # Summarise the copy number of each AR122 and BAC120 markers.
        tln_summary_file = TlnTableSummaryFile(outdir, prefix)
        ar122_copy_number_file = CopyNumberFileAR122(outdir, prefix)
        bac120_copy_number_file = CopyNumberFileBAC120(outdir, prefix)

        # Process each genome.
        for db_genome_id, info in tqdm_log(sorted(gene_dict.items()), unit='genome'):
            cur_marker_dir = os.path.join(outdir, DIR_MARKER_GENE)
            pfam_tophit_file = TopHitPfamFile(cur_marker_dir, db_genome_id)
            tigr_tophit_file = TopHitTigrFile(cur_marker_dir, db_genome_id)
            pfam_tophit_file.read()
            tigr_tophit_file.read()

            # Summarise each of the markers for this genome.
            ar122_copy_number_file.add_genome(db_genome_id, info.get("aa_gene_path"),
                                              pfam_tophit_file, tigr_tophit_file)
            bac120_copy_number_file.add_genome(db_genome_id, info.get("aa_gene_path"),
                                               pfam_tophit_file, tigr_tophit_file)

            # Write the best translation table to disk for this genome.
            tln_summary_file.add_genome(db_genome_id, info.get("best_translation_table"))

        # Write each of the summary files to disk.
        ar122_copy_number_file.write()
        bac120_copy_number_file.write()
        tln_summary_file.write()

        # Create a symlink to store the summary files in the root.
        symlink_f(PATH_BAC120_MARKER_SUMMARY.format(prefix=prefix),
                  os.path.join(outdir, os.path.basename(PATH_BAC120_MARKER_SUMMARY.format(prefix=prefix))))
        symlink_f(PATH_AR122_MARKER_SUMMARY.format(prefix=prefix),
                  os.path.join(outdir, os.path.basename(PATH_AR122_MARKER_SUMMARY.format(prefix=prefix))))
        symlink_f(PATH_TLN_TABLE_SUMMARY.format(prefix=prefix),
                  os.path.join(outdir, os.path.basename(PATH_TLN_TABLE_SUMMARY.format(prefix=prefix))))

        # Write the single copy AR122/BAC120 FASTA files to disk.
        if write_single_copy_genes:
            fasta_dir = os.path.join(outdir, DIR_IDENTIFY_FASTA)
            self.logger.info(f'Writing unaligned single-copy genes to: {fasta_dir}')

            # Iterate over each domain.
            marker_doms = list()
            marker_doms.append((Config.AR122_MARKERS['PFAM'] +
                                Config.AR122_MARKERS['TIGRFAM'],
                                ar122_copy_number_file, 'ar122'))
            marker_doms.append((Config.BAC120_MARKERS['PFAM'] +
                                Config.BAC120_MARKERS['TIGRFAM'],
                                bac120_copy_number_file, 'bac120'))
            for marker_names, marker_file, marker_d in marker_doms:

                # Create the domain-specific subdirectory.
                fasta_d_dir = os.path.join(fasta_dir, marker_d)
                make_sure_path_exists(fasta_d_dir)

                # Iterate over each marker.
                for marker_name in marker_names:
                    marker_name = marker_name.rstrip(r'\.[HMMhmm]')
                    marker_path = os.path.join(fasta_d_dir, f'{marker_name}.fa')

                    to_write = list()
                    for genome_id in sorted(gene_dict):
                        unq_hits = marker_file.get_single_copy_hits(genome_id)
                        if marker_name in unq_hits:
                            to_write.append(f'>{genome_id}')
                            to_write.append(unq_hits[marker_name]['seq'])

                    if len(to_write) > 0:
                        with open(marker_path, 'w') as fh:
                            fh.write('\n'.join(to_write))

    def identify(self, genomes, tln_tables, out_dir, prefix, force, write_single_copy_genes):
        """Identify marker genes in genomes.

        Parameters
        ----------
        genomes : dict
            Genome IDs as the key, path to genome file as value.
        tln_tables: Dict[str, int]
            Genome ID -> translation table mapping for those user-specified.
        out_dir : str
            Path to the output directory.
        prefix : str
            Prefix to append to generated files.
        force : bool
            Overwrite any existing files.
        write_single_copy_genes : bool
            Write unique AR122/BAC120 marker files to disk.

        Raises
        ------
        GTDBTkException
            If an exception is encountered during the identify step.

        """
        check_dependencies(['prodigal', 'hmmsearch'])

        self.logger.info(f'Identifying markers in {len(genomes):,} genomes with '
                         f'{self.cpus} threads.')

        self.marker_gene_dir = os.path.join(out_dir, DIR_MARKER_GENE)
        prodigal = Prodigal(self.cpus,
                            self.marker_gene_dir,
                            self.protein_file_suffix,
                            self.nt_gene_file_suffix,
                            self.gff_file_suffix,
                            force)
        self.logger.log(Config.LOG_TASK, f'Running Prodigal {prodigal.version} to identify genes.')
        genome_dictionary = prodigal.run(genomes, tln_tables)

        # annotated genes against TIGRFAM and Pfam databases
        self.logger.log(Config.LOG_TASK, 'Identifying TIGRFAM protein families.')
        gene_files = [genome_dictionary[db_genome_id]['aa_gene_path']
                      for db_genome_id in genome_dictionary.keys()]
        tigr_search = TigrfamSearch(self.cpus,
                                    self.tigrfam_hmms,
                                    self.protein_file_suffix,
                                    self.tigrfam_suffix,
                                    self.tigrfam_top_hit_suffix,
                                    self.checksum_suffix,
                                    self.marker_gene_dir)
        tigr_search.run(gene_files)

        self.logger.log(Config.LOG_TASK, 'Identifying Pfam protein families.')
        pfam_search = PfamSearch(self.cpus,
                                 self.pfam_hmm_dir,
                                 self.protein_file_suffix,
                                 self.pfam_suffix,
                                 self.pfam_top_hit_suffix,
                                 self.checksum_suffix,
                                 self.marker_gene_dir)
        pfam_search.run(gene_files)
        self.logger.info(f'Annotations done using HMMER {tigr_search.version}.')

        self.logger.log(Config.LOG_TASK, 'Summarising identified marker genes.')
        self._report_identified_marker_genes(genome_dictionary, out_dir, prefix,
                                             write_single_copy_genes)

    def _path_to_identify_data(self, identity_dir, warn=True):
        """Get path to genome data produced by 'identify' command."""

        marker_gene_dir = os.path.join(identity_dir, DIR_MARKER_GENE)

        genomic_files = {}
        lq_gids = list()
        for gid in os.listdir(marker_gene_dir):
            gid_dir = os.path.join(marker_gene_dir, gid)
            if not os.path.isdir(gid_dir):
                continue

            aa_gene_path = os.path.join(
                gid_dir, gid + self.protein_file_suffix)

            # Check if any genes were called
            if os.path.getsize(aa_gene_path) < 1:
                lq_gids.append(gid)
            else:
                genomic_files[gid] = {'aa_gene_path': aa_gene_path,
                                      'translation_table_path': TlnTableFile.get_path(gid_dir, gid),
                                      'nt_gene_path': os.path.join(gid_dir, gid + self.nt_gene_file_suffix),
                                      'gff_path': os.path.join(gid_dir, gid + self.gff_file_suffix)
                                      }

        if len(lq_gids) > 0 and warn:
            self.logger.warning(f'Excluding {len(lq_gids)} genomes '
                                f'in the identify directory which have no genes '
                                f'called (see gtdb.warnings.log)')
            self.warnings.warning(f'Excluding the following {len(lq_gids)} genomes '
                                  f'which were found in the identify directory '
                                  f'with no genes called.')
            for lq_gid in lq_gids:
                self.warnings.info(lq_gid)
        return genomic_files

    def _msa_filter_by_taxa(self, concatenated_file: str,
                            gtdb_taxonomy: Dict[str, Tuple[str, str, str, str, str, str, str]],
                            taxa_filter: Optional[str],
                            outgroup_taxon: Optional[str]) -> Dict[str, str]:
        """Filter GTDB MSA to a subset of specified taxa.

        Parameters
        ----------
        concatenated_file
            The path to the MSA.
        gtdb_taxonomy
            A dictionary mapping the accession to the 7 rank taxonomy.
        taxa_filter
            A comma separated list of taxa to include.
        outgroup_taxon
            If using an outgroup (de novo workflow), ensure this is retained.

        Returns
        -------
        Dict[str, str]
            The genome id to msa of those genomes specified in the filter.
        """

        msa = read_fasta(concatenated_file)
        msa_len = len(msa)
        self.logger.info(f'Read concatenated alignment for {msa_len:,} GTDB genomes.')

        if taxa_filter is not None:
            taxa_to_keep = set(taxa_filter.split(','))

            if outgroup_taxon not in taxa_to_keep and outgroup_taxon is not None:
                taxa_to_keep.add(outgroup_taxon)

            filtered_genomes = 0
            for genome_id, taxa in gtdb_taxonomy.items():
                common_taxa = taxa_to_keep.intersection(taxa)
                if len(common_taxa) == 0:
                    if genome_id in msa:
                        del msa[genome_id]
                        filtered_genomes += 1

            msg = f'Filtered {filtered_genomes / msa_len:.2%} ({filtered_genomes:,}/{msa_len:,}) ' \
                  f'taxa based on assigned taxonomy, {msa_len - filtered_genomes:,} taxa remain.'
            self.logger.info(msg) if len(msa) > 0 else self.logger.warning(msg)

        return msa

    def _apply_mask(self, gtdb_msa, user_msa, msa_mask, min_perc_aa):
        """Apply canonical mask to MSA file."""
        aligned_genomes = merge_two_dicts(gtdb_msa, user_msa)
        list_mask = np.fromfile(msa_mask, dtype='S1') == b'1'

        output_seqs, pruned_seqs = dict(), dict()
        for seq_id, seq in tqdm_log(aligned_genomes.items(), unit='sequence'):
            list_seq = np.fromiter(seq, dtype='S1')
            if list_mask.shape[0] != list_seq.shape[0]:
                raise MSAMaskLengthMismatch(
                    'Mask and alignment length do not match.')

            list_masked_seq = list_seq[list_mask]

            masked_seq_unique = np.unique(list_masked_seq, return_counts=True)
            masked_seq_counts = defaultdict(lambda: 0)
            for aa_char, aa_count in zip(masked_seq_unique[0], masked_seq_unique[1]):
                masked_seq_counts[aa_char.decode('utf-8')] = aa_count

            masked_seq = list_masked_seq.tostring().decode('utf-8')

            valid_bases = list_masked_seq.shape[0] - masked_seq_counts['.'] - masked_seq_counts['-']
            if seq_id in user_msa and valid_bases < list_masked_seq.shape[0] * min_perc_aa:
                pruned_seqs[seq_id] = masked_seq
                continue

            output_seqs[seq_id] = masked_seq

        return output_seqs, pruned_seqs

    def _write_msa(self, seqs, output_file, gtdb_taxonomy):
        """Write sequences to FASTA file."""

        with open(output_file, 'w') as fout:
            for genome_id, alignment in sorted(seqs.items()):
                if genome_id in gtdb_taxonomy:
                    fout.write('>%s %s\n' %
                               (genome_id, ';'.join(gtdb_taxonomy[genome_id])))
                else:
                    fout.write('>%s\n' % genome_id)
                fout.write('%s\n' % alignment)

    def genome_domain(self, identity_dir, prefix):
        """Determine domain of User genomes based on identified marker genes."""
        bac_count = defaultdict(int)
        ar_count = defaultdict(int)

        # Load the marker files for each domain
        ar122_marker_file = CopyNumberFileAR122(identity_dir, prefix)
        ar122_marker_file.read()
        bac120_marker_file = CopyNumberFileBAC120(identity_dir, prefix)
        bac120_marker_file.read()

        # Get the number of single copy markers for each domain
        for out_d, marker_summary in ((ar_count, ar122_marker_file),
                                      (bac_count, bac120_marker_file)):
            for genome_id in marker_summary.genomes:
                out_d[genome_id] = len(marker_summary.get_single_copy_hits(genome_id))

        bac_gids = set()
        ar_gids = set()
        bac_ar_diff = {}
        for gid in bac_count:
            arc_aa_per = (ar_count[gid] * 100.0 / Config.AR_MARKER_COUNT)
            bac_aa_per = (bac_count[gid] * 100.0 / Config.BAC_MARKER_COUNT)
            if bac_aa_per >= arc_aa_per:
                bac_gids.add(gid)
            else:
                ar_gids.add(gid)
            if abs(bac_aa_per - arc_aa_per) <= 10:
                bac_ar_diff[gid] = {'bac120': round(
                    bac_aa_per, 1), 'ar122': round(arc_aa_per, 1)}

        return bac_gids, ar_gids, bac_ar_diff

    def _write_marker_info(self, marker_db, marker_file):
        """Write out information about markers comprising MSA."""

        marker_paths = {"PFAM": os.path.join(self.pfam_hmm_dir, 'individual_hmms'),
                        "TIGRFAM": os.path.join(os.path.dirname(self.tigrfam_hmms), 'individual_hmms')}

        with open(marker_file, 'w') as fh:
            fh.write('Marker ID\tName\tDescription\tLength (bp)\n')
            for db_marker in sorted(marker_db):
                for marker in marker_db[db_marker]:
                    marker_id = marker[0:marker.rfind('.')]
                    marker_path = os.path.join(marker_paths[db_marker], marker)

                    # get marker name, description, and size
                    with open(marker_path) as fp:
                        for line in fp:
                            if line.startswith("NAME  "):
                                marker_name = line.split("  ")[1].strip()
                            elif line.startswith("DESC  "):
                                marker_desc = line.split("  ")[1].strip()
                            elif line.startswith("LENG  "):
                                marker_size = line.split("  ")[1].strip()
                                break
                    fh.write('%s\t%s\t%s\t%s\n' %
                             (marker_id, marker_name, marker_desc, marker_size))

    def align(self,
              identify_dir,
              skip_gtdb_refs,
              taxa_filter,
              min_perc_aa,
              custom_msa_filters,
              skip_trimming,
              rnd_seed,
              cols_per_gene,
              min_consensus,
              max_consensus,
              min_per_taxa,
              out_dir,
              prefix,
              outgroup_taxon,
              genomes_to_process=None):
        """Align marker genes in genomes."""

        # If the user is re-running this step, check if the identify step is consistent.
        genomic_files = self._path_to_identify_data(identify_dir, identify_dir != out_dir)
        if genomes_to_process is not None and len(genomic_files) != len(genomes_to_process):
            self.logger.error('{} are not present in the input list of genome to process.'.format(
                list(set(genomic_files.keys()) - set(genomes_to_process.keys()))))
            raise InconsistentGenomeBatch(
                'You are attempting to run GTDB-Tk on a non-empty directory that contains extra '
                'genomes not present in your initial identify directory. Remove them, or run '
                'GTDB-Tk on a new directory.')

        # If this is being run as a part of classify_wf, copy the required files.
        if identify_dir != out_dir:
            identify_path = os.path.join(out_dir, DIR_IDENTIFY)
            make_sure_path_exists(identify_path)
            copy(CopyNumberFileBAC120(identify_dir, prefix).path, identify_path)
            copy(CopyNumberFileAR122(identify_dir, prefix).path, identify_path)
            copy(TlnTableSummaryFile(identify_dir, prefix).path, identify_path)

        # Create the align intermediate directory.
        make_sure_path_exists(os.path.join(out_dir, DIR_ALIGN_INTERMEDIATE))

        # Write out files with marker information
        ar122_marker_info_file = MarkerInfoFileAR122(out_dir, prefix)
        ar122_marker_info_file.write()
        bac120_marker_info_file = MarkerInfoFileBAC120(out_dir, prefix)
        bac120_marker_info_file.write()

        # Determine what domain each genome belongs to.
        bac_gids, ar_gids, _bac_ar_diff = self.genome_domain(identify_dir, prefix)

        # # Create a temporary directory that will be used to generate each of the alignments.
        # with tempfile.TemporaryDirectory(prefix='gtdbtk_tmp_') as dir_tmp_arc, \
        #         tempfile.TemporaryDirectory(prefix='gtdbtk_tmp_') as dir_tmp_bac:
        #
        #     cur_gid_dict = {x: genomic_files[x] for x in ar_gids}
        #     self.logger.info(f'Collecting marker sequences from {len(cur_gid_dict):,} '
        #                      f'genomes identified as archaeal.')
        #     align.concat_single_copy_hits(dir_tmp_arc,
        #                                   cur_gid_dict,
        #                                   ar122_marker_info_file)
        #

        self.logger.info(f'Aligning markers in {len(genomic_files):,} genomes with {self.cpus} CPUs.')
        dom_iter = ((bac_gids, Config.CONCAT_BAC120, Config.MASK_BAC120, "bac120", 'bacterial', CopyNumberFileBAC120),
                    (ar_gids, Config.CONCAT_AR122, Config.MASK_AR122, "ar122", 'archaeal', CopyNumberFileAR122))
        gtdb_taxonomy = Taxonomy().read(self.taxonomy_file)
        for gids, msa_file, mask_file, marker_set_id, domain_str, copy_number_f in dom_iter:

            # No genomes identified as this domain.
            if len(gids) == 0:
                continue

            self.logger.info(f'Processing {len(gids):,} genomes identified as {domain_str}.')
            if marker_set_id == 'bac120':
                marker_info_file = bac120_marker_info_file
                marker_filtered_genomes = os.path.join(
                    out_dir, PATH_BAC120_FILTERED_GENOMES.format(prefix=prefix))
                marker_msa_path = os.path.join(
                    out_dir, PATH_BAC120_MSA.format(prefix=prefix))
                marker_user_msa_path = os.path.join(
                    out_dir, PATH_BAC120_USER_MSA.format(prefix=prefix))
            else:
                marker_info_file = ar122_marker_info_file
                marker_filtered_genomes = os.path.join(
                    out_dir, PATH_AR122_FILTERED_GENOMES.format(prefix=prefix))
                marker_msa_path = os.path.join(
                    out_dir, PATH_AR122_MSA.format(prefix=prefix))
                marker_user_msa_path = os.path.join(
                    out_dir, PATH_AR122_USER_MSA.format(prefix=prefix))

            cur_genome_files = {
                gid: f for gid, f in genomic_files.items() if gid in gids}

            if skip_gtdb_refs:
                gtdb_msa = {}
            else:
                gtdb_msa = self._msa_filter_by_taxa(msa_file,
                                                    gtdb_taxonomy,
                                                    taxa_filter,
                                                    outgroup_taxon)
            gtdb_msa_mask = os.path.join(Config.MASK_DIR, mask_file)

            # Generate the user MSA.
            user_msa = align.align_marker_set(cur_genome_files, marker_info_file, copy_number_f, self.cpus)

            # self.logger.log(Config.LOG_TASK, f'Aligning {len(cur_genome_files):,} {domain_str} genomes.')
            # hmm_aligner = HmmAligner(self.cpus,
            #                          self.pfam_top_hit_suffix,
            #                          self.tigrfam_top_hit_suffix,
            #                          self.protein_file_suffix,
            #                          self.pfam_hmm_dir,
            #                          self.tigrfam_hmms,
            #                          Config.BAC120_MARKERS,
            #                          Config.AR122_MARKERS)
            # user_msa = hmm_aligner.align_marker_set(cur_genome_files,
            #                                         marker_set_id)

            # Write the individual marker alignments to disk
            if self.debug:
                self._write_individual_markers(
                    user_msa, marker_set_id, marker_info_file.path, out_dir, prefix)

            # filter columns without sufficient representation across taxa
            if skip_trimming:
                self.logger.info('Skipping custom filtering and selection of columns.')
                pruned_seqs = {}
                trimmed_seqs = merge_two_dicts(gtdb_msa, user_msa)

            elif custom_msa_filters:
                aligned_genomes = merge_two_dicts(gtdb_msa, user_msa)
                self.logger.info('Performing custom filtering and selection of columns.')

                trim_msa = TrimMSA(cols_per_gene,
                                   min_perc_aa / 100.0,
                                   min_consensus / 100.0,
                                   max_consensus / 100.0,
                                   min_per_taxa / 100.0,
                                   rnd_seed,
                                   os.path.join(out_dir, f'filter_{marker_set_id}'))

                trimmed_seqs, pruned_seqs = trim_msa.trim(aligned_genomes,
                                                          marker_info_file.path)

                if trimmed_seqs:
                    self.logger.info('Filtered MSA from {:,} to {:,} AAs.'.format(
                        len(list(aligned_genomes.values())[0]),
                        len(list(trimmed_seqs.values())[0])))

                self.logger.info('Filtered {:,} genomes with amino acids in <{:.1f}% of columns in filtered MSA.'.format(
                    len(pruned_seqs),
                    min_perc_aa))

                filtered_user_genomes = set(
                    pruned_seqs).intersection(user_msa)
                if len(filtered_user_genomes):
                    self.logger.info('Filtered genomes include {:.} user submitted genomes.'.format(len(
                        filtered_user_genomes)))
            else:
                self.logger.log(Config.LOG_TASK,
                    f'Masking columns of {domain_str} multiple sequence alignment using canonical mask.')
                trimmed_seqs, pruned_seqs = self._apply_mask(gtdb_msa,
                                                             user_msa,
                                                             gtdb_msa_mask,
                                                             min_perc_aa / 100.0)
                self.logger.info('Masked {} alignment from {:,} to {:,} AAs.'.format(
                                    domain_str,
                                    len(list(user_msa.values())[0]),
                                    len(list(trimmed_seqs.values())[0])))

                if min_perc_aa > 0:
                    self.logger.info('{:,} {} user genomes have amino acids in <{:.1f}% of columns in filtered MSA.'.format(
                        len(pruned_seqs),
                        domain_str,
                        min_perc_aa))

            # write out filtering information
            with open(marker_filtered_genomes, 'w') as fout:
                for pruned_seq_id, pruned_seq in pruned_seqs.items():
                    if len(pruned_seq) == 0:
                        perc_alignment = 0
                    else:
                        valid_bases = sum([1 for c in pruned_seq if c.isalpha()])
                        perc_alignment = valid_bases * 100.0 / len(pruned_seq)
                    fout.write(f'{pruned_seq_id}\tInsufficient number of amino acids in MSA ({perc_alignment:.1f}%)\n')

            # write out MSAs
            if not skip_gtdb_refs:
                self.logger.info(f'Creating concatenated alignment for {len(trimmed_seqs):,} '
                                 f'{domain_str} GTDB and user genomes.')
                self._write_msa(trimmed_seqs, marker_msa_path, gtdb_taxonomy)

            trimmed_user_msa = {k: v for k, v in trimmed_seqs.items()
                                if k in user_msa}
            if len(trimmed_user_msa) > 0:
                self.logger.info(f'Creating concatenated alignment for {len(trimmed_user_msa):,} '
                                 f'{domain_str} user genomes.')
                self._write_msa(trimmed_user_msa,
                                marker_user_msa_path, gtdb_taxonomy)
            else:
                self.logger.info(f'All {domain_str} user genomes have been filtered out.')

            # Create symlinks to the summary files
            if marker_set_id == 'bac120':
                symlink_f(PATH_BAC120_FILTERED_GENOMES.format(prefix=prefix),
                          os.path.join(out_dir, os.path.basename(PATH_BAC120_FILTERED_GENOMES.format(prefix=prefix))))
                if len(trimmed_user_msa) > 0:
                    symlink_f(PATH_BAC120_USER_MSA.format(prefix=prefix),
                              os.path.join(out_dir, os.path.basename(PATH_BAC120_USER_MSA.format(prefix=prefix))))
                if not skip_gtdb_refs:
                    symlink_f(PATH_BAC120_MSA.format(prefix=prefix),
                              os.path.join(out_dir, os.path.basename(PATH_BAC120_MSA.format(prefix=prefix))))
            elif marker_set_id == 'ar122':
                symlink_f(PATH_AR122_FILTERED_GENOMES.format(prefix=prefix),
                          os.path.join(out_dir, os.path.basename(PATH_AR122_FILTERED_GENOMES.format(prefix=prefix))))
                if len(trimmed_user_msa) > 0:
                    symlink_f(PATH_AR122_USER_MSA.format(prefix=prefix),
                              os.path.join(out_dir, os.path.basename(PATH_AR122_USER_MSA.format(prefix=prefix))))
                if not skip_gtdb_refs:
                    symlink_f(PATH_AR122_MSA.format(prefix=prefix),
                              os.path.join(out_dir, os.path.basename(PATH_AR122_MSA.format(prefix=prefix))))
            else:
                raise GenomeMarkerSetUnknown('There was an error determining the marker set.')

    def _write_individual_markers(self, user_msa, marker_set_id, marker_list, out_dir, prefix):
        marker_dir = join(out_dir, DIR_ALIGN_MARKERS)
        make_sure_path_exists(marker_dir)

        markers, total_msa_len = self._parse_marker_info_file(marker_list)
        marker_to_msa = dict()
        offset = 0
        for marker_id, marker_desc, marker_len in sorted(markers, key=lambda x: x[0]):
            path_msa = os.path.join(marker_dir, f'{prefix}.{marker_set_id}.{marker_id}.faa')
            marker_to_msa[path_msa] = defaultdict(str)
            for gid, msa in user_msa.items():
                marker_to_msa[path_msa][gid] += msa[offset: marker_len + offset]
            offset += marker_len

        if total_msa_len != offset:
            self.logger.warning(
                'Internal error: the total MSA length is not equal to the offset.')
        for path_marker, gid_dict in marker_to_msa.items():
            with open(path_marker, 'w') as fh:
                for genome_id, genome_msa in gid_dict.items():
                    fh.write(f'>{genome_id}\n{genome_msa}\n')
        self.logger.debug(f'Successfully written all markers to: {marker_dir}')

    def _parse_marker_info_file(self, path):
        total_msa_len = 0
        markers = list()
        with open(path, 'r') as fh:
            fh.readline()
            for line in fh.readlines():
                list_info = line.split("\t")
                marker_id = list_info[0]
                marker_name = '%s: %s' % (list_info[1], list_info[2])
                marker_len = int(list_info[3])
                markers.append((marker_id, marker_name, marker_len))
                total_msa_len += marker_len
        return markers, total_msa_len
