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
from collections import defaultdict
from shutil import copy

import numpy as np

import config.config as Config
from biolib_lite.execute import check_dependencies
from biolib_lite.seq_io import read_fasta
from biolib_lite.taxonomy import Taxonomy
from external.hmm_aligner import HmmAligner
from external.pfam_search import PfamSearch
from external.prodigal import Prodigal
from external.tigrfam_search import TigrfamSearch
from gtdbtk.config.output import *
from gtdbtk.exceptions import GenomeMarkerSetUnknown, MSAMaskLengthMismatch
from tools import merge_two_dicts, symlink_f
from trim_msa import TrimMSA


class Markers(object):
    """Identify and align marker genes."""

    def __init__(self, cpus=1):
        """Initialize."""

        self.logger = logging.getLogger('timestamp')

        self.cpus = cpus

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

    def _report_identified_marker_genes(self, gene_dict, outdir, marker_gene_dir, prefix):
        """Report statistics for identified marker genes."""

        translation_table_file = open(os.path.join(outdir, PATH_TLN_TABLE_SUMMARY.format(prefix=prefix)), "w")
        bac_outfile = open(os.path.join(outdir, PATH_BAC120_MARKER_SUMMARY.format(prefix=prefix)), "w")
        arc_outfile = open(os.path.join(outdir, PATH_AR122_MARKER_SUMMARY.format(prefix=prefix)), "w")

        header = "Name\tnumber_unique_genes\tnumber_multiple_genes\tnumber_missing_genes\tlist_unique_genes\tlist_multiple_genes\tlist_missing_genes\n"

        bac_outfile.write(header)
        arc_outfile.write(header)

        # gather information for all marker genes
        marker_dbs = {"PFAM": PFAM_TOP_HIT_SUFFIX, "TIGR": TIGRFAM_TOP_HIT_SUFFIX}

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
                tophit_path = protein_file.replace(PROTEIN_FILE_SUFFIX, marker_suffix)

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

                            if (markerid not in marker_bac_list_original and
                                    markerid not in marker_arc_list_original):
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

            translation_table_file.write('{}\t{}\n'.format(
                db_genome_id, info.get("best_translation_table")))

        bac_outfile.close()
        arc_outfile.close()
        translation_table_file.close()

        # Create a symlink to store the summary files in the root.
        symlink_f(PATH_BAC120_MARKER_SUMMARY.format(prefix=prefix),
                  os.path.join(outdir, os.path.basename(PATH_BAC120_MARKER_SUMMARY.format(prefix=prefix))))
        symlink_f(PATH_AR122_MARKER_SUMMARY.format(prefix=prefix),
                  os.path.join(outdir, os.path.basename(PATH_AR122_MARKER_SUMMARY.format(prefix=prefix))))
        symlink_f(PATH_TLN_TABLE_SUMMARY.format(prefix=prefix),
                  os.path.join(outdir, os.path.basename(PATH_TLN_TABLE_SUMMARY.format(prefix=prefix))))

    def identify(self, genomes, out_dir, prefix, force):
        """Identify marker genes in genomes.

        Parameters
        ----------
        genomes : dict
            Genome IDs as the key, path to genome file as value.
        out_dir : str
            Path to the output directory.
        prefix : str
            Prefix to append to generated files.
        force : bool
            Overwrite any existing files.

        Raises
        ------
        GTDBTkException
            If an exception is encountered during the identify step.

        """
        check_dependencies(['prodigal', 'hmmsearch'])

        self.logger.info('Identifying markers in %d genomes with %d threads.' % (len(genomes),
                                                                                 self.cpus))

        self.logger.info("Running Prodigal to identify genes.")
        self.marker_gene_dir = os.path.join(out_dir, DIR_MARKER_GENE)
        prodigal = Prodigal(self.cpus,
                            False,
                            self.marker_gene_dir,
                            self.protein_file_suffix,
                            self.nt_gene_file_suffix,
                            self.gff_file_suffix,
                            force)
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
            genome_dictionary, out_dir, self.marker_gene_dir, prefix)

    def _path_to_identify_data(self, identity_dir):
        """Get path to genome data produced by 'identify' command."""

        marker_gene_dir = os.path.join(identity_dir, DIR_MARKER_GENE)

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
        msa_len = len(msa)
        self.logger.info(
            'Read concatenated alignment for %d GTDB genomes.' % msa_len)

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

            msg = 'Filtered %.2f%% (%d/%d) taxa based on assigned taxonomy, ' \
                  '%d taxa remain.' % (
                      (float(filtered_genomes) / float(msa_len)) * 100.0,
                      filtered_genomes, msa_len, msa_len - filtered_genomes)
            self.logger.info(msg) if len(msa) > 0 else self.logger.warning(msg)

        return msa

    def _apply_mask(self, gtdb_msa, user_msa, msa_mask, min_perc_aa):
        """Apply canonical mask to MSA file."""

        aligned_genomes = merge_two_dicts(gtdb_msa, user_msa)

        with open(msa_mask, 'r') as f:
            mask = f.readline().strip()
        list_mask = np.array([True if c == '1' else False for c in mask], dtype=bool)

        if len(mask) != len(aligned_genomes.values()[0]):
            self.logger.error('Mask and alignment length do not match.')
            raise MSAMaskLengthMismatch

        output_seqs = {}
        pruned_seqs = {}
        for seq_id, seq in aligned_genomes.iteritems():
            masked_seq = ''.join(np.array(list(seq), dtype=str)[list_mask])

            valid_bases = len(masked_seq) - masked_seq.count('.') - masked_seq.count('-')
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

    def genome_domain(self, identity_dir, prefix):
        """Determine domain of User genomes based on identified marker genes."""

        bac_count = defaultdict(int)
        ar_count = defaultdict(int)

        for d, marker_file in ((bac_count, os.path.join(identity_dir, PATH_BAC120_MARKER_SUMMARY.format(prefix=prefix))),
                               (ar_count, os.path.join(identity_dir, PATH_AR122_MARKER_SUMMARY.format(prefix=prefix)))):
            with open(marker_file) as f:
                f.readline()

                for line in f:
                    line_split = line.strip().split('\t')

                    gid = line_split[0]
                    num_markers = int(line_split[1])

                    d[gid] = num_markers

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

        fout = open(marker_file, 'w')
        fout.write('Marker ID\tName\tDescription\tLength (bp)\n')
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

                fout.write('%s\t%s\t%s\t%s\n' %
                           (marker_id, marker_name, marker_desc, marker_size))
        fout.close()

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
              outgroup_taxon):
        """Align marker genes in genomes."""

        if identify_dir != out_dir:
            if not os.path.isdir(os.path.join(out_dir, DIR_IDENTIFY)):
                os.makedirs(os.path.join(out_dir, DIR_IDENTIFY))

            copy(os.path.join(identify_dir, PATH_BAC120_MARKER_SUMMARY.format(prefix=prefix)),
                 os.path.join(out_dir, DIR_IDENTIFY))
            copy(os.path.join(identify_dir, PATH_AR122_MARKER_SUMMARY.format(prefix=prefix)),
                 os.path.join(out_dir, DIR_IDENTIFY))

            identify_gene_file = os.path.join(identify_dir, PATH_TLN_TABLE_SUMMARY.format(prefix=prefix))
            copy(identify_gene_file, os.path.join(out_dir, DIR_IDENTIFY))

        if not os.path.exists(os.path.join(out_dir, DIR_ALIGN_INTERMEDIATE)):
            os.makedirs(os.path.join(out_dir, DIR_ALIGN_INTERMEDIATE))

        # write out files with marker information
        bac120_marker_info_file = os.path.join(out_dir, PATH_BAC120_MARKER_INFO.format(prefix=prefix))
        self._write_marker_info(Config.BAC120_MARKERS, bac120_marker_info_file)
        ar122_marker_info_file = os.path.join(out_dir, PATH_AR122_MARKER_INFO.format(prefix=prefix))
        self._write_marker_info(Config.AR122_MARKERS, ar122_marker_info_file)

        genomic_files = self._path_to_identify_data(identify_dir)
        self.logger.info('Aligning markers in %d genomes with %d threads.' % (len(genomic_files),
                                                                              self.cpus))

        # determine marker set for each user genome
        bac_gids, ar_gids, _bac_ar_diff = self.genome_domain(identify_dir, prefix)

        # align user genomes
        gtdb_taxonomy = Taxonomy().read(self.taxonomy_file)
        for gids, msa_file, mask_file, marker_set_id in ((bac_gids, Config.CONCAT_BAC120, Config.MASK_BAC120, "bac120"),
                                                         (ar_gids, Config.CONCAT_AR122, Config.MASK_AR122, "ar122")):

            if len(gids) == 0:
                continue

            if marker_set_id == 'bac120':
                self.logger.info('Processing %d genomes identified as bacterial.' % len(gids))
                marker_info_file = bac120_marker_info_file
                marker_filtered_genomes = os.path.join(out_dir, PATH_BAC120_FILTERED_GENOMES.format(prefix=prefix))
                marker_msa_path = os.path.join(out_dir, PATH_BAC120_MSA.format(prefix=prefix))
                marker_user_msa_path = os.path.join(out_dir, PATH_BAC120_USER_MSA.format(prefix=prefix))
            else:
                self.logger.info('Processing %d genomes identified as archaeal.' % len(gids))
                marker_info_file = ar122_marker_info_file
                marker_filtered_genomes = os.path.join(out_dir, PATH_AR122_FILTERED_GENOMES.format(prefix=prefix))
                marker_msa_path = os.path.join(out_dir, PATH_AR122_MSA.format(prefix=prefix))
                marker_user_msa_path = os.path.join(out_dir, PATH_AR122_USER_MSA.format(prefix=prefix))

            cur_genome_files = {
                gid: f for gid, f in genomic_files.iteritems() if gid in gids}

            if skip_gtdb_refs:
                gtdb_msa = {}
            else:
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
            if skip_trimming:
                self.logger.info(
                    'Skipping custom filtering and selection of columns.')
                pruned_seqs = {}
                trimmed_seqs = merge_two_dicts(gtdb_msa, user_msa)

            elif custom_msa_filters:
                aligned_genomes = merge_two_dicts(gtdb_msa, user_msa)
                self.logger.info(
                    'Performing custom filtering and selection of columns.')

                trim_msa = TrimMSA(cols_per_gene,
                                   min_perc_aa / 100.0,
                                   min_consensus / 100.0,
                                   max_consensus / 100.0,
                                   min_per_taxa / 100.0,
                                   rnd_seed,
                                   os.path.join(out_dir, 'filter_%s' % marker_set_id))

                trimmed_seqs, pruned_seqs = trim_msa.trim(aligned_genomes,
                                                          marker_info_file)

                if trimmed_seqs:
                    self.logger.info('Filtered MSA from %d to %d AAs.' % (
                        len(aligned_genomes.values()[0]),
                        len(trimmed_seqs.values()[0])))

                self.logger.info('Filtered %d genomes with amino acids in <%.1f%% of columns in filtered MSA.' % (
                    len(pruned_seqs),
                    min_perc_aa))

                filtered_user_genomes = set(
                    pruned_seqs).intersection(user_msa)
                if len(filtered_user_genomes):
                    self.logger.info('Filtered genomes include %d user submitted genomes.' % len(
                        filtered_user_genomes))
            else:
                self.logger.info(
                    'Masking columns of multiple sequence alignment using canonical mask.')
                trimmed_seqs, pruned_seqs = self._apply_mask(gtdb_msa,
                                                             user_msa,
                                                             gtdb_msa_mask,
                                                             min_perc_aa / 100.0)
                self.logger.info('Masked alignment from %d to %d AAs.' % (len(user_msa.values()[0]),
                                                                          len(trimmed_seqs.values()[0])))

                if min_perc_aa > 0:
                    self.logger.info('%d user genomes have amino acids in <%.1f%% of columns in filtered MSA.' % (
                        len(pruned_seqs),
                        min_perc_aa))

            # write out filtering information
            with open(marker_filtered_genomes, 'w') as fout:
                for pruned_seq_id, pruned_seq in pruned_seqs.items():
                    if len(pruned_seq) == 0:
                        perc_alignment = 0
                    else:
                        valid_bases = sum(
                            [1 for c in pruned_seq if c.isalpha()])
                        perc_alignment = valid_bases * 100.0 / len(pruned_seq)
                    fout.write('%s\t%s\n' % (pruned_seq_id,
                                             'Insufficient number of amino acids in MSA (%.1f%%)' % perc_alignment))

            # write out MSAs
            if not skip_gtdb_refs:
                self.logger.info('Creating concatenated alignment for %d GTDB and user genomes.' % len(trimmed_seqs))
                self._write_msa(trimmed_seqs, marker_msa_path, gtdb_taxonomy)

            trimmed_user_msa = {
                k: v for k, v in trimmed_seqs.iteritems() if k in user_msa}
            if len(trimmed_user_msa) > 0:
                self.logger.info('Creating concatenated alignment for %d user genomes.' % len(trimmed_user_msa))
                self._write_msa(trimmed_user_msa, marker_user_msa_path, gtdb_taxonomy)
            else:
                if marker_set_id == 'bac120':
                    self.logger.info('All bacterial user genomes have been filtered out.')
                else:
                    self.logger.info('All archaeal user genomes have been filtered out.')

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
                self.logger.error('There was an error determining the marker set.')
                raise GenomeMarkerSetUnknown
