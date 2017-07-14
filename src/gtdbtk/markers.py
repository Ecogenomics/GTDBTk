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
import tempfile
import shutil
import multiprocessing
import time
import subprocess

from gtdblib.trimming import trim_seqs

from biolib.common import remove_extension
from biolib.seq_io import read_fasta
from biolib.taxonomy import Taxonomy

from external.prodigal import Prodigal
from external.tigrfam_search import TigrfamSearch
from external.pfam_search import PfamSearch
from external.hmm_aligner import HmmAligner

import config.config as Config
import config.config_metadata as ConfigMetadata

from tools import splitchunks, list_genomes_dir, merge_two_dicts, genomes_to_process


class Markers(object):
    """Identify and align marker genes."""

    def __init__(self, cpus=1):
        """Initialize."""
        
        self.logger = logging.getLogger('timestamp')
        
        self.cpus = cpus
        
        self.identify_dir = 'genome_data'

        self.genome_file_suffix = ConfigMetadata.GENOME_FILE_SUFFIX
        self.marker_gene_dir = Config.MARKER_GENE_DIR
        self.protein_file_suffix = ConfigMetadata.PROTEIN_FILE_SUFFIX
        self.nt_gene_file_suffix = ConfigMetadata.NT_GENE_FILE_SUFFIX
        self.gff_file_suffix = ConfigMetadata.GFF_FILE_SUFFIX
        self.checksum_suffix = ConfigMetadata.CHECKSUM_SUFFIX
        
        self.taxonomy_file = Config.TAXONOMY_FILE

        self.pfam_hmm_dir = ConfigMetadata.PFAM_HMM_DIR
        self.pfam_suffix = ConfigMetadata.PFAM_SUFFIX
        self.pfam_top_hit_suffix = ConfigMetadata.PFAM_TOP_HIT_SUFFIX

        self.tigrfam_hmms = ConfigMetadata.TIGRFAM_HMMS
        self.tigrfam_suffix = ConfigMetadata.TIGRFAM_SUFFIX
        self.tigrfam_top_hit_suffix = ConfigMetadata.TIGRFAM_TOP_HIT_SUFFIX
                
    def _prepare_genomes(self, dict_user_genomes, output_dir):
        """Copy genome files in temporary folders in order to process them.

        Parameters
        ----------
        dict_user_genomes : str
            Name of file describing genomes to add.
        output_dir : str
            Output directory.

        Returns
        -------
        dict
            Dictionary indicating the temporary gene file for each genome.
        """

        genomic_files = {}
        for user_genome_raw in dict_user_genomes.keys():
            user_genome = os.path.splitext(os.path.basename(user_genome_raw))[0]
            genome_output_dir = os.path.join(output_dir, user_genome)
            if not os.path.exists(genome_output_dir):
                os.makedirs(genome_output_dir)
                
            prodigal_output_dir = os.path.join(genome_output_dir, self.marker_gene_dir)
            if not os.path.exists(prodigal_output_dir):
                os.makedirs(prodigal_output_dir)

            # copy the genome file to the temporary folder
            fasta_target_file = os.path.join(genome_output_dir, user_genome + self.genome_file_suffix)
            shutil.copy(dict_user_genomes.get(user_genome_raw), fasta_target_file)

            genomic_files[user_genome] = fasta_target_file
            
        return genomic_files
        
    def _report_identified_marker_genes(self, gene_dict, outdir, prefix):
        """Report statistics for identified marker genes."""

        bac_outfile = open(os.path.join(outdir, prefix + "_bac120_markers_summary.tsv"), "w")
        arc_outfile = open(os.path.join(outdir, prefix + "_ar122_markers_summary.tsv"), "w")
        rps23_outfile = open(os.path.join(outdir, prefix + "_rps23_markers_summary.tsv"), "w")

        header = "Name\tnumber_unique_genes\tnumber_multiple_genes\tnumber_missing_genes\tlist_unique_genes\tlist_multiple_genes\tlist_missing_genes\n"

        bac_outfile.write(header)
        arc_outfile.write(header)
        rps23_outfile.write(header)

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
         
        marker_rps23_list_original = []
        for db_marker in Config.RPS23_MARKERS.keys():
            marker_rps23_list_original.extend([marker.replace(".HMM", "").replace(".hmm", "") 
                                                for marker in Config.RPS23_MARKERS[db_marker]])

        gene_bac_list, gene_arc_list, gene_rps23_list = [], [], []

        for db_genome_id, info in gene_dict.items():
            unique_genes_bac, multi_hits_bac, missing_genes_bac = [], [], []
            unique_genes_arc, multi_hits_arc, missing_genes_arc = [], [], []
            unique_genes_rps23, multi_hits_rps23, missing_genes_rps23 = [], [], []

            gene_bac_dict, gene_arc_dict, gene_rps23_dict = {}, {}, {}

            path = info.get("aa_gene_path")
            for marker_db, marker_suffix in marker_dbs.iteritems():
                # get all gene sequences
                protein_file = str(path)
                tophit_path = protein_file.replace(ConfigMetadata.PROTEIN_FILE_SUFFIX, marker_suffix)
                
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
                                and markerid not in marker_arc_list_original
                                and markerid not in marker_rps23_list_original):
                                continue

                            if markerid in marker_bac_list_original:
                                if markerid in gene_bac_dict:
                                    gene_bac_dict.get(markerid)["multihit"] = True
                                else:
                                    gene_bac_dict[markerid] = {
                                        "gene": genename,
                                        #"gene_seq": all_genes_dict.get(genename),
                                        "multihit": False}

                            if markerid in marker_arc_list_original:
                                if markerid in gene_arc_dict:
                                    gene_arc_dict.get(markerid)["multihit"] = True
                                else:
                                    gene_arc_dict[markerid] = {
                                        "gene": genename,
                                        #"gene_seq": all_genes_dict.get(genename),
                                        "multihit": False}
                                        
                            if markerid in marker_rps23_list_original:
                                if markerid in gene_rps23_dict:
                                    gene_rps23_dict.get(markerid)["multihit"] = True
                                else:
                                    gene_rps23_dict[markerid] = {
                                        "gene": genename,
                                        #"gene_seq": all_genes_dict.get(genename),
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
                    
            for mid in marker_rps23_list_original:
                if mid not in gene_rps23_dict:
                    missing_genes_rps23.append(mid)
                elif gene_rps23_dict[mid]["multihit"]:
                    multi_hits_rps23.append(mid)
                else:
                    unique_genes_rps23.append(mid)

            bac_outfile.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format(db_genome_id,
                                                                           len(unique_genes_bac),
                                                                           len(multi_hits_bac),
                                                                           len(missing_genes_bac),
                                                                           ','.join(unique_genes_bac),
                                                                           ','.join(multi_hits_bac),
                                                                           ','.join(missing_genes_bac)))

            arc_outfile.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format(db_genome_id,
                                                                           len(unique_genes_arc),
                                                                           len(multi_hits_arc),
                                                                           len(missing_genes_arc),
                                                                           ','.join(unique_genes_arc),
                                                                           ','.join(multi_hits_arc),
                                                                           ','.join(missing_genes_arc)))
                                                                           
            rps23_outfile.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format(db_genome_id,
                                                                           len(unique_genes_rps23),
                                                                           len(multi_hits_rps23),
                                                                           len(missing_genes_rps23),
                                                                           ','.join(unique_genes_rps23),
                                                                           ','.join(multi_hits_rps23),
                                                                           ','.join(missing_genes_rps23)))

        bac_outfile.close()
        arc_outfile.close()
        rps23_outfile.close()

    def identify(self, genome_dir, batchfile, proteins, out_dir, prefix):
        """Identify marker genes in genomes."""
        
        try:
            genomes = genomes_to_process(genome_dir, batchfile)
            self.logger.info('Identifying markers in %d genomes with %d threads.' % (len(genomes), 
                                                                                        self.cpus))
            
            # gather information for all marker genes
            output_dir = os.path.join(out_dir, self.identify_dir)
            genomic_files = self._prepare_genomes(genomes, 
                                                    output_dir)
 
            self.logger.info("Running Prodigal to identify genes.")
            prodigal = Prodigal(self.cpus,
                                proteins,
                                self.marker_gene_dir, 
                                self.protein_file_suffix,
                                self.nt_gene_file_suffix, 
                                self.gff_file_suffix)
            genome_dictionary = prodigal.run(genomic_files)

            # annotated genes against TIGRfam and Pfam databases
            self.logger.info("Identifying TIGRfam protein families.")
            gene_files = [genome_dictionary[db_genome_id]['aa_gene_path']
                          for db_genome_id in genome_dictionary.keys()]

            tigr_search = TigrfamSearch(self.cpus, 
                                        self.tigrfam_hmms, 
                                        self.protein_file_suffix,
                                        self.tigrfam_suffix, 
                                        self.tigrfam_top_hit_suffix, 
                                        self.checksum_suffix)
            tigr_search.run(gene_files)

            self.logger.info("Identifying Pfam protein families.")
            pfam_search = PfamSearch(self.cpus, 
                                        self.pfam_hmm_dir, 
                                        self.protein_file_suffix,
                                        self.pfam_suffix, 
                                        self.pfam_top_hit_suffix, 
                                        self.checksum_suffix)
            pfam_search.run(gene_files)

            self._report_identified_marker_genes(genome_dictionary, out_dir, prefix)

        except IOError as e:
            self.logger.error(str(e))
            self.logger.error("GTDB-Tk has encountered an error.")
            
        except Exception as e:
            self.logger.error(str(e))
            raise
            
    def _path_to_identify_data(self, genome_ids, indir):
        """Get path to genome data produced by 'identify' command."""
                
        genomic_files = {}
        for genome_id_raw in genome_ids:
            genome_id = os.path.splitext(os.path.basename(genome_id_raw))[0]
            genomic_files[genome_id] = {'aa_gene_path': os.path.join(indir, self.identify_dir, genome_id, 'marker_genes', genome_id + self.protein_file_suffix),
                                   'translation_table_path': os.path.join(indir, self.identify_dir, genome_id, 'marker_genes', 'prodigal_translation_table.tsv'),
                                   'nt_gene_path': os.path.join(indir, self.identify_dir, genome_id, 'marker_genes', genome_id + self.nt_gene_file_suffix),
                                   'gff_path': os.path.join(indir, self.identify_dir, genome_id, 'marker_genes', genome_id + self.gff_file_suffix)
                                   }
        return genomic_files
        
    def _msa_filter_by_taxa(self, concatenated_file, gtdb_taxonomy, taxa_filter):
        """Filter GTDB MSA filtered to specified taxa."""
        
        msa = read_fasta(concatenated_file)
        self.logger.info('Read concatenated alignment for %d GTDB genomes.' % len(msa))
        
        if taxa_filter is not None:
            taxa_to_keep = set(taxa_filter.split(','))

            filtered_genomes = 0
            for genome_id, taxa in gtdb_taxonomy.iteritems():
                common_taxa = taxa_to_keep.intersection(taxa)
                if len(common_taxa) == 0:
                    if genome_id in msa:
                        del msa[genome_id]
                        filtered_genomes += 1
                        
            self.logger.info('Filtered %d taxa based on assigned taxonomy.' % filtered_genomes)
        
        return msa
        
    def _apply_mask(self, alignments, msa_mask, min_perc_aa):
        """Apply canonical mask to MSA file."""
        
        mask = open(msa_mask).readline().strip()
        
        if len(mask) != len(alignments.values()[0]):
            self.logger.error('Mask and alignment length do not match.')
            sys.exit()
        
        output_seqs = {}
        pruned_seqs = {}
        for seq_id, seq in alignments.iteritems():
            masked_seq = ''.join([seq[i] for i in xrange(0, len(mask)) if mask[i] == '1'])
            
            valid_bases = len(masked_seq) - masked_seq.count('.') - masked_seq.count('-')
            if valid_bases < len(masked_seq) * min_perc_aa:
                pruned_seqs[seq_id] = masked_seq
                continue

            output_seqs[seq_id] = masked_seq

        return output_seqs, pruned_seqs
        
    def _write_msa(self, seqs, output_file, gtdb_taxonomy):
        """Write sequences to FASTA file."""
        
        fout = open(output_file, 'w')
        for genome_id, alignment in seqs.iteritems():
            if genome_id in gtdb_taxonomy:
                fout.write('>%s %s\n' % (genome_id, ';'.join(gtdb_taxonomy[genome_id])))
            else:
                fout.write('>%s\n' % genome_id)
            fout.write('%s\n' % alignment)
        fout.close()

    def align(self,
                genome_dir,
                batchfile, 
                identify_dir, 
                marker_set_id, 
                taxa_filter, 
                min_perc_aa,
                custom_msa_filters,
                consensus, 
                min_per_taxa, 
                out_dir, 
                prefix):
        """Align marker genes in genomes."""

        try:
            genomes = genomes_to_process(genome_dir, batchfile)
            self.logger.info('Aligning markers in %d genomes with %d threads.' % (len(genomes), 
                                                                                        self.cpus))
            genomic_files = self._path_to_identify_data(genomes, identify_dir)
            
            gtdb_taxonomy = Taxonomy().read(self.taxonomy_file)
                                                                                
            if marker_set_id == 'bac120':
                gtdb_msa = self._msa_filter_by_taxa(Config.CONCAT_BAC120,
                                                        gtdb_taxonomy,
                                                        taxa_filter)
                gtdb_msa_mask = os.path.join(Config.MASK_DIR, Config.MASK_BAC120)
            elif marker_set_id == 'ar122':
                gtdb_msa = self._msa_filter_by_taxa(Config.CONCAT_AR122,
                                                        gtdb_taxonomy,
                                                        taxa_filter)
                gtdb_msa_mask = os.path.join(Config.MASK_DIR, Config.MASK_AR122)
            elif marker_set_id == 'rps23':
                gtdb_msa = self._msa_filter_by_taxa(Config.CONCAT_RPS23,
                                                        gtdb_taxonomy,
                                                        taxa_filter)
                gtdb_msa_mask = os.path.join(Config.MASK_DIR, Config.MASK_RPS23)
            else:
                self.logger.error('Unrecognized marker set: %s' % marker_set_id)
                sys.exit()

            hmm_aligner = HmmAligner(self.cpus,
                                        self.pfam_top_hit_suffix,
                                        self.tigrfam_top_hit_suffix,
                                        self.protein_file_suffix,
                                        self.pfam_hmm_dir,
                                        self.tigrfam_hmms,
                                        Config.BAC120_MARKERS,
                                        Config.AR122_MARKERS,
                                        Config.RPS23_MARKERS)
            user_msa = hmm_aligner.align_marker_set(genomic_files, 
                                                    marker_set_id)

            # filter columns without sufficient representation across taxa
            aligned_genomes = merge_two_dicts(gtdb_msa, user_msa)
            if custom_msa_filters:
                self.logger.info('Trimming columns with insufficient taxa or poor consensus.')
                trimmed_seqs, pruned_seqs, count_wrong_pa, count_wrong_cons = trim_seqs(aligned_genomes, 
                                                                                        min_per_taxa / 100.0, 
                                                                                        consensus / 100.0, 
                                                                                        min_perc_aa / 100.0)
                self.logger.info(('Trimmed alignment from %d to %d AA (%d by minimum taxa percent, '
                                    + '%d by consensus).') % (len(aligned_genomes.values()[0]),
                                                                len(trimmed_seqs.values()[0]), 
                                                                count_wrong_pa, 
                                                                count_wrong_cons))
            else:
                self.logger.info('Masking columns of multiple sequence alignment.')
                trimmed_seqs, pruned_seqs = self._apply_mask(aligned_genomes, 
                                                                gtdb_msa_mask, 
                                                                min_perc_aa / 100.0)
                self.logger.info('Masked alignment from %d to %d AA.' % (len(aligned_genomes.values()[0]),
                                                                            len(trimmed_seqs.values()[0])))                                                                
                                                            
            self.logger.info('Pruned %d taxa with amino acids in <%.1f%% of columns in filtered MSA.' % (
                                len(pruned_seqs), 
                                min_perc_aa))
                                
            pruned_user_genomes = set(pruned_seqs).intersection(user_msa)
            if len(pruned_user_genomes):
                self.logger.info('Pruned genomes include %d user submitted genomes.' % len(pruned_user_genomes))

            # write out MSA
            self.logger.info('Creating concatenated alignment for %d taxa.' % len(trimmed_seqs)) 
            msa_file = os.path.join(out_dir, prefix + ".msa.faa")
            self._write_msa(trimmed_seqs, msa_file, gtdb_taxonomy)
            
            user_msa_file = os.path.join(out_dir, prefix + ".user_msa.faa")
            trimmed_user_msa = {k:v for k, v in trimmed_seqs.iteritems() if k in user_msa}
            self._write_msa(trimmed_user_msa, user_msa_file, gtdb_taxonomy)

        except IOError as e:
            self.logger.error(str(e))
            self.logger.error("GTDB-Tk has encountered an error.")
