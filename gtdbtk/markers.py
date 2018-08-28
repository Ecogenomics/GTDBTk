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

import os
import sys
import logging
from collections import defaultdict

from biolib.seq_tk import trim_seqs

from biolib.common import remove_extension
from biolib.seq_io import read_fasta
from biolib.taxonomy import Taxonomy
from biolib.external.execute import check_dependencies

from external.prodigal import Prodigal
from external.tigrfam_search import TigrfamSearch
from external.pfam_search import PfamSearch
from external.hmm_aligner import HmmAligner

import config.config as Config
import config.config_metadata as ConfigMetadata

from tools import merge_two_dicts


class Markers(object):
    """Identify and align marker genes."""

    def __init__(self, cpus=1):
        """Initialize."""

        self.logger = logging.getLogger('timestamp')

        self.cpus = cpus

        self.genome_file_suffix = ConfigMetadata.GENOME_FILE_SUFFIX
        self.protein_file_suffix = ConfigMetadata.PROTEIN_FILE_SUFFIX
        self.nt_gene_file_suffix = ConfigMetadata.NT_GENE_FILE_SUFFIX
        self.gff_file_suffix = ConfigMetadata.GFF_FILE_SUFFIX
        self.checksum_suffix = ConfigMetadata.CHECKSUM_SUFFIX

        self.taxonomy_file = Config.TAXONOMY_FILE

        self.pfam_hmm_dir = Config.PFAM_HMM_DIR
        self.pfam_suffix = ConfigMetadata.PFAM_SUFFIX
        self.pfam_top_hit_suffix = ConfigMetadata.PFAM_TOP_HIT_SUFFIX

        self.tigrfam_hmms = Config.TIGRFAM_HMMS
        self.tigrfam_suffix = ConfigMetadata.TIGRFAM_SUFFIX
        self.tigrfam_top_hit_suffix = ConfigMetadata.TIGRFAM_TOP_HIT_SUFFIX

    def _report_identified_marker_genes(self, gene_dict, outdir, prefix):
        """Report statistics for identified marker genes."""

        bac_outfile = open(os.path.join(
            outdir, prefix + "_bac120_markers_summary.tsv"), "w")
        arc_outfile = open(os.path.join(
            outdir, prefix + "_ar122_markers_summary.tsv"), "w")

        header = "Name\tnumber_unique_genes\tnumber_multiple_genes\tnumber_missing_genes\tlist_unique_genes\tlist_multiple_genes\tlist_missing_genes\n"

        bac_outfile.write(header)
        arc_outfile.write(header)

        # gather information for all marker genes
        marker_dbs = {"PFAM": ConfigMetadata.PFAM_TOP_HIT_SUFFIX,
                      "TIGR": ConfigMetadata.TIGRFAM_TOP_HIT_SUFFIX}

        marker_bac_list_original = []
        for db_marker in Config.BAC120_MARKERS.keys():
            marker_bac_list_original.extend([marker.replace(".HMM", "").replace(".hmm", "")
                                             for marker in Config.BAC120_MARKERS[db_marker]])

        marker_arc_list_original = []
        for db_marker in Config.AR122_MARKERS.keys():
            marker_arc_list_original.extend([marker.replace(".HMM", "").replace(".hmm", "")
                                             for marker in Config.AR122_MARKERS[db_marker]])

        for db_genome_id, info in gene_dict.items():
            unique_genes_bac, multi_hits_bac, missing_genes_bac = [], [], []
            unique_genes_arc, multi_hits_arc, missing_genes_arc = [], [], []

            gene_bac_dict, gene_arc_dict = {}, {}

            path = info.get("aa_gene_path")
            for _marker_db, marker_suffix in marker_dbs.iteritems():
                # get all gene sequences
                protein_file = str(path)
                tophit_path = protein_file.replace(
                    ConfigMetadata.PROTEIN_FILE_SUFFIX, marker_suffix)

                # we load the list of all the genes detected in the genome
                all_genes_dict = read_fasta(protein_file, False)

                # Prodigal adds an asterisks at the end of each called genes.
                # These asterisks sometimes appear in the MSA, which can be
                # an issue for some downstream software
                for seq_id, seq in all_genes_dict.iteritems():
                    if seq[-1] == '*':
                        all_genes_dict[seq_id] = seq[:-1]

                # we store the tophit file line by line and store the
                # information in a dictionary
                with open(tophit_path) as tp:
                    # first line is header line
                    tp.readline()

                    for line_tp in tp:
                        linelist = line_tp.split("\t")
                        genename = linelist[0]
                        sublist = linelist[1]
                        if ";" in sublist:
                            diff_markers = sublist.split(";")
                        else:
                            diff_markers = [sublist]

                        for each_mark in diff_markers:
                            sublist = each_mark.split(",")
                            markerid = sublist[0]

                            if (markerid not in marker_bac_list_original
                                    and markerid not in marker_arc_list_original):
                                continue

                            if markerid in marker_bac_list_original:
                                if markerid in gene_bac_dict:
                                    gene_bac_dict.get(markerid)[
                                        "multihit"] = True
                                else:
                                    gene_bac_dict[markerid] = {
                                        "gene": genename,
                                        "multihit": False}

                            if markerid in marker_arc_list_original:
                                if markerid in gene_arc_dict:
                                    gene_arc_dict.get(markerid)[
                                        "multihit"] = True
                                else:
                                    gene_arc_dict[markerid] = {
                                        "gene": genename,
                                        "multihit": False}

            for mid in marker_bac_list_original:
                if mid not in gene_bac_dict:
                    missing_genes_bac.append(mid)
                elif gene_bac_dict[mid]["multihit"]:
                    multi_hits_bac.append(mid)
                else:
                    unique_genes_bac.append(mid)

            for mid in marker_arc_list_original:
                if mid not in gene_arc_dict:
                    missing_genes_arc.append(mid)
                elif gene_arc_dict[mid]["multihit"]:
                    multi_hits_arc.append(mid)
                else:
                    unique_genes_arc.append(mid)

            bac_outfile.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format(db_genome_id,
                                                                           len(unique_genes_bac),
                                                                           len(multi_hits_bac),
                                                                           len(missing_genes_bac),
                                                                           ','.join(
                                                                               unique_genes_bac),
                                                                           ','.join(
                                                                               multi_hits_bac),
                                                                           ','.join(missing_genes_bac)))

            arc_outfile.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format(db_genome_id,
                                                                           len(unique_genes_arc),
                                                                           len(multi_hits_arc),
                                                                           len(missing_genes_arc),
                                                                           ','.join(
                                                                               unique_genes_arc),
                                                                           ','.join(
                                                                               multi_hits_arc),
                                                                           ','.join(missing_genes_arc)))

        bac_outfile.close()
        arc_outfile.close()

    def identify(self,
                 genomes,
                 out_dir,
                 prefix):
        """Identify marker genes in genomes."""

        check_dependencies(['prodigal', 'hmmsearch'])

        try:
            self.logger.info('Identifying markers in %d genomes with %d threads.' % (len(genomes),
                                                                                     self.cpus))

            self.logger.info("Running Prodigal to identify genes.")
            self.marker_gene_dir = os.path.join(
                out_dir, Config.MARKER_GENE_DIR)
            prodigal = Prodigal(self.cpus,
                                False,
                                self.marker_gene_dir,
                                self.protein_file_suffix,
                                self.nt_gene_file_suffix,
                                self.gff_file_suffix)
            genome_dictionary = prodigal.run(genomes)

            # annotated genes against TIGRFAM and Pfam databases
            self.logger.info("Identifying TIGRFAM protein families.")
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

            self.logger.info("Identifying Pfam protein families.")
            pfam_search = PfamSearch(self.cpus,
                                     self.pfam_hmm_dir,
                                     self.protein_file_suffix,
                                     self.pfam_suffix,
                                     self.pfam_top_hit_suffix,
                                     self.checksum_suffix,
                                     self.marker_gene_dir)
            pfam_search.run(gene_files)

            self._report_identified_marker_genes(
                genome_dictionary, out_dir, prefix)

        except IOError as e:
            self.logger.error(str(e))
            self.logger.error("GTDB-Tk has encountered an error.")

        except Exception as e:
            self.logger.error(str(e))
            raise

    def _path_to_identify_data(self, identity_dir):
        """Get path to genome data produced by 'identify' command."""

        marker_gene_dir = os.path.join(identity_dir, Config.MARKER_GENE_DIR)

        genomic_files = {}
        for gid in os.listdir(marker_gene_dir):
            gid_dir = os.path.join(marker_gene_dir, gid)
            if not os.path.isdir(gid_dir):
                continue

            genomic_files[gid] = {'aa_gene_path': os.path.join(gid_dir, gid + self.protein_file_suffix),
                                  'translation_table_path': os.path.join(gid_dir, 'prodigal_translation_table.tsv'),
                                  'nt_gene_path': os.path.join(gid_dir, gid + self.nt_gene_file_suffix),
                                  'gff_path': os.path.join(gid_dir, gid + self.gff_file_suffix)
                                  }
        return genomic_files

    def _msa_filter_by_taxa(self, concatenated_file, gtdb_taxonomy, taxa_filter, outgroup_taxon):
        """Filter GTDB MSA filtered to specified taxa."""

        msa = read_fasta(concatenated_file)
        self.logger.info(
            'Read concatenated alignment for %d GTDB genomes.' % len(msa))

        if taxa_filter is not None:
            taxa_to_keep = set(taxa_filter.split(','))

            if outgroup_taxon not in taxa_to_keep and outgroup_taxon is not None:
                taxa_to_keep.add(outgroup_taxon)

            filtered_genomes = 0
            for genome_id, taxa in gtdb_taxonomy.iteritems():
                common_taxa = taxa_to_keep.intersection(taxa)
                if len(common_taxa) == 0:
                    if genome_id in msa:
                        del msa[genome_id]
                        filtered_genomes += 1

            self.logger.info(
                'Filtered %d taxa based on assigned taxonomy.' % filtered_genomes)

        return msa

    def _apply_mask(self, gtdb_msa, user_msa, msa_mask, min_perc_aa):
        """Apply canonical mask to MSA file."""

        aligned_genomes = merge_two_dicts(gtdb_msa, user_msa)

        mask = open(msa_mask).readline().strip()

        if len(mask) != len(aligned_genomes.values()[0]):
            self.logger.error('Mask and alignment length do not match.')
            sys.exit()

        output_seqs = {}
        pruned_seqs = {}
        for seq_id, seq in aligned_genomes.iteritems():
            masked_seq = ''.join(
                [seq[i] for i in xrange(0, len(mask)) if mask[i] == '1'])

            valid_bases = len(masked_seq) - \
                masked_seq.count('.') - masked_seq.count('-')
            if seq_id in user_msa and valid_bases < len(masked_seq) * min_perc_aa:
                pruned_seqs[seq_id] = masked_seq
                continue

            output_seqs[seq_id] = masked_seq

        return output_seqs, pruned_seqs

    def _write_msa(self, seqs, output_file, gtdb_taxonomy):
        """Write sequences to FASTA file."""

        fout = open(output_file, 'w')
        for genome_id, alignment in seqs.iteritems():
            if genome_id in gtdb_taxonomy:
                fout.write('>%s %s\n' %
                           (genome_id, ';'.join(gtdb_taxonomy[genome_id])))
            else:
                fout.write('>%s\n' % genome_id)
            fout.write('%s\n' % alignment)
        fout.close()

    def _genome_domain(self, identity_dir, prefix):
        """Determine domain of User genomes based on identified marker genes."""

        bac_count = defaultdict(int)
        ar_count = defaultdict(int)
        for d, marker_file in ((bac_count, prefix + '_bac120_markers_summary.tsv'),
                               (ar_count, prefix + '_ar122_markers_summary.tsv')):
            with open(os.path.join(identity_dir, marker_file)) as f:
                f.readline()

                for line in f:
                    line_split = line.strip().split('\t')

                    gid = line_split[0]
                    num_markers = int(line_split[1])

                    d[gid] = num_markers

        bac_gids = set()
        ar_gids = set()
        for gid in bac_count:
            arc_aa_per = (ar_count[gid] * 100.0 / Config.AR_MARKER_COUNT)
            bac_aa_per = (bac_count[gid] * 100.0 / Config.BAC_MARKER_COUNT)
            if bac_aa_per >= arc_aa_per:
                bac_gids.add(gid)
            else:
                ar_gids.add(gid)

        return bac_gids, ar_gids

    def align(self,
              identify_dir,
              taxa_filter,
              min_perc_aa,
              custom_msa_filters,
              consensus,
              min_per_taxa,
              out_dir,
              prefix,
              outgroup_taxon):
        """Align marker genes in genomes."""

        try:
            genomic_files = self._path_to_identify_data(identify_dir)
            self.logger.info('Aligning markers in %d genomes with %d threads.' % (len(genomic_files),
                                                                                  self.cpus))

            # determine marker set for each user genome
            bac_gids, ar_gids = self._genome_domain(identify_dir, prefix)

            # align user genomes
            gtdb_taxonomy = Taxonomy().read(self.taxonomy_file)
            for gids, msa_file, mask_file, marker_set_id in ((bac_gids, Config.CONCAT_BAC120, Config.MASK_BAC120, "bac120"),
                                                             (ar_gids, Config.CONCAT_AR122, Config.MASK_AR122, "ar122")):

                if len(gids) == 0:
                    continue

                if marker_set_id == 'bac120':
                    self.logger.info(
                        'Processing %d genomes identified as bacterial.' % len(gids))
                else:
                    self.logger.info(
                        'Processing %d genomes identified as archaeal.' % len(gids))

                cur_genome_files = {
                    gid: f for gid, f in genomic_files.iteritems() if gid in gids}

                gtdb_msa = self._msa_filter_by_taxa(msa_file,
                                                    gtdb_taxonomy,
                                                    taxa_filter,
                                                    outgroup_taxon)
                gtdb_msa_mask = os.path.join(Config.MASK_DIR, mask_file)

                hmm_aligner = HmmAligner(self.cpus,
                                         self.pfam_top_hit_suffix,
                                         self.tigrfam_top_hit_suffix,
                                         self.protein_file_suffix,
                                         self.pfam_hmm_dir,
                                         self.tigrfam_hmms,
                                         Config.BAC120_MARKERS,
                                         Config.AR122_MARKERS,
                                         Config.RPS23_MARKERS)
                user_msa = hmm_aligner.align_marker_set(cur_genome_files,
                                                        marker_set_id)

                # filter columns without sufficient representation across taxa
                if custom_msa_filters:
                    aligned_genomes = merge_two_dicts(gtdb_msa, user_msa)
                    self.logger.info(
                        'Trimming columns with insufficient taxa or poor consensus.')
                    trimmed_seqs, pruned_seqs, count_wrong_pa, count_wrong_cons = trim_seqs(aligned_genomes,
                                                                                            min_per_taxa / 100.0,
                                                                                            consensus / 100.0,
                                                                                            min_perc_aa / 100.0)
                    self.logger.info(('Trimmed alignment from %d to %d AA (%d by minimum taxa percent, ' +
                                      '%d by consensus).') % (len(aligned_genomes.values()[0]),
                                                              len(trimmed_seqs.values()[
                                                                  0]),
                                                              count_wrong_pa,
                                                              count_wrong_cons))

                    self.logger.info('%d taxa have amino acids in <%.1f%% of columns in filtered MSA.' % (
                        len(pruned_seqs),
                        min_perc_aa))

                    pruned_user_genomes = set(
                        pruned_seqs).intersection(user_msa)
                    if len(pruned_user_genomes):
                        self.logger.info(
                            'Pruned genomes include %d user submitted genomes.' % len(pruned_user_genomes))
                else:
                    self.logger.info(
                        'Masking columns of multiple sequence alignment.')
                    trimmed_seqs, pruned_seqs = self._apply_mask(gtdb_msa,
                                                                 user_msa,
                                                                 gtdb_msa_mask,
                                                                 min_perc_aa / 100.0)
                    self.logger.info('Masked alignment from %d to %d AA.' % (len(gtdb_msa.values()[0]),
                                                                             len(trimmed_seqs.values()[0])))

                    if min_perc_aa > 0:
                        self.logger.info('%d user genomes have amino acids in <%.1f%% of columns in filtered MSA.' % (
                            len(pruned_seqs),
                            min_perc_aa))

                # write out filtering information
                fout = open(os.path.join(out_dir, prefix +
                                         ".%s.filtered.tsv" % marker_set_id), 'w')
                for pruned_seq_id, pruned_seq in pruned_seqs.items():
                    valid_bases = len(
                        pruned_seq) - pruned_seq.count('.') - pruned_seq.count('-')
                    perc_alignment = valid_bases * 100.0 / len(pruned_seq)
                    fout.write('%s\t%s\n' % (pruned_seq_id,
                                             'Insufficient number of amino acids in MSA (%.1f%%)' % perc_alignment))
                fout.close()

                # write out MSA
                self.logger.info(
                    'Creating concatenated alignment for %d taxa.' % len(trimmed_seqs))
                msa_file = os.path.join(
                    out_dir, prefix + ".%s.msa.fasta" % marker_set_id)
                self._write_msa(trimmed_seqs, msa_file, gtdb_taxonomy)

                user_msa_file = os.path.join(
                    out_dir, prefix + ".%s.user_msa.fasta" % marker_set_id)
                trimmed_user_msa = {
                    k: v for k, v in trimmed_seqs.iteritems() if k in user_msa}
                self._write_msa(trimmed_user_msa, user_msa_file, gtdb_taxonomy)

                #==============================================================
                # #all_user_msa_file = os.path.join(out_dir, prefix + ".%s.user_msa.fasta" % marker_set_id)
                # trimmed_all_user_msa = {k:v for k, v in trimmed_seqs.iteritems() if k in user_msa}
                # pruned_all_user_msa = {k:v for k, v in pruned_seqs.iteritems() if k in user_msa}
                # all_user_msa = merge_two_dicts(trimmed_all_user_msa,pruned_all_user_msa)
                # self._write_msa(all_user_msa, all_user_msa_file, gtdb_taxonomy)
                #==============================================================

        except IOError as e:
            self.logger.error(str(e))
            self.logger.error("GTDB-Tk has encountered an error.")
