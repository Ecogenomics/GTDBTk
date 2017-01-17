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


from biolib.seq_io import read_fasta
from gtdblib.trimming import trim_seqs
from biolib.external.fasttree import FastTree

from external.Prodigal import Prodigal
from external.TigrfamSearch import TigrfamSearch
from external.PfamSearch import PfamSearch

import config.Config as Config
import config.ConfigMetadata as ConfigMetadata
import config.DefaultValues as DefaultValues

from Tools import fastaPathGenerator, splitchunks, list_genomes_dir, merge_two_dicts


class GtdbManager(object):
    """Manages genomes, concatenated alignment, and metadata for tree inference and visualization."""

    def __init__(self, threads=1):
        self.init = True
        self.threads = threads

        self.genomeFileSuffix = ConfigMetadata.GENOME_FILE_SUFFIX
        self.userAnnotationDir = Config.ANNOTATION_DIR
        self.protein_file_suffix = ConfigMetadata.PROTEIN_FILE_SUFFIX
        self.nt_gene_file_suffix = ConfigMetadata.NT_GENE_FILE_SUFFIX
        self.gff_file_suffix = ConfigMetadata.GFF_FILE_SUFFIX
        self.checksum_suffix = ConfigMetadata.CHECKSUM_SUFFIX

        self.concatenated_bacteria = Config.CONCAT_BAC
        self.concatenated_archaea = Config.CONCAT_ARC

        self.taxonomy_file = Config.TAXONOMY_FILE

        self.pfam_hmm_dir = ConfigMetadata.PFAM_HMM_DIR
        self.pfam_suffix = ConfigMetadata.PFAM_SUFFIX
        self.pfam_top_hit_suffix = ConfigMetadata.PFAM_TOP_HIT_SUFFIX

        self.tigrfam_hmms = ConfigMetadata.TIGRFAM_HMMS
        self.tigrfam_suffix = ConfigMetadata.TIGRFAM_SUFFIX
        self.tigrfam_top_hit_suffix = ConfigMetadata.TIGRFAM_TOP_HIT_SUFFIX

        self.logger = logging.getLogger()
        self.logger.setLevel(logging.DEBUG)

        ch = logging.StreamHandler(sys.stdout)
        ch.setLevel(logging.DEBUG)
        formatter = logging.Formatter('%(asctime)s - <%(levelname)s>: %(message)s')
        ch.setFormatter(formatter)
        self.logger.addHandler(ch)

    def MakeTreeData(self, batchfile, tree_generate, outfile=None, maindomain=None, filter_taxa=None, min_perc_aa=None,
                     consensus=None, min_per_taxa=None, outdir=None, notree=None, prefix=None):
        """Generate a Tree based on the MSA.

        Parameters
        ----------
        user_directory : str
            Directory path where all fasta files are.
        tree_generate : boolean
            if True : a tree is created. if False, the domain of each genome is predicted
        outfile : str
            Path to the output file for domain prediction
        maindomain: str
            "bacteria" or "archaea". This flag informa which marker set and which MSA to use.
        filter_taxa: str
            Select the taxa of interest to generate the tree
        """

        # create tree data
        try:

            self.tmp_output_dir = tempfile.mkdtemp()
            print self.tmp_output_dir

            dict_user_genomes = self.list_genomes(batchfile)
            genomic_files = self._prepareGenomes(dict_user_genomes, self.tmp_output_dir)

            self.logger.info("Running Prodigal to identify genes.")
            prodigal = Prodigal(self.threads, self.userAnnotationDir, self.protein_file_suffix,
                                self.nt_gene_file_suffix, self.gff_file_suffix, self.checksum_suffix)
            genome_dictionary = prodigal.run(genomic_files)

            # annotated genes against TIGRfam and Pfam databases
            self.logger.info("Identifying TIGRfam protein families.")
            gene_files = [genome_dictionary[db_genome_id]['aa_gene_path']
                          for db_genome_id in genome_dictionary.keys()]
            tigr_search = TigrfamSearch(self.threads, self.tigrfam_hmms, self.protein_file_suffix,
                                        self.tigrfam_suffix, self.tigrfam_top_hit_suffix, self.checksum_suffix)
            tigr_search.run(gene_files)

            self.logger.info("Identifying Pfam protein families.")
            pfam_search = PfamSearch(self.threads, self.pfam_hmm_dir, self.protein_file_suffix,
                                     self.pfam_suffix, self.pfam_top_hit_suffix, self.checksum_suffix)
            pfam_search.run(gene_files)

            taxonomy_dict = {}
            if maindomain == 'bacteria':
                dict_gtdb_genomes = self._selectGTDBGenomes(self.concatenated_bacteria, filter_taxa)
            else:
                dict_gtdb_genomes = self._selectGTDBGenomes(self.concatenated_archaea, filter_taxa)

            # if tree_generate:
            # we generate a tree based on the MSA
            # if not:
            # we do not generate the tree and only predict the domain of the potential genome based on the frequencies of markers.
            if tree_generate:
                dict_aligned_genomes = self._alignMarkerSet(genome_dictionary, maindomain)
            else:
                dict_aligned_genomes = self._calculatePresenceMarkerSets(genome_dictionary, outfile)
                return True

            # filter columns without sufficient representation across taxa
            self.logger.info('Trimming columns with insufficient taxa or poor consensus.')

            dict_aligned_genomes = merge_two_dicts(dict_gtdb_genomes, dict_aligned_genomes)
            trimmed_seqs, pruned_seqs, count_wrong_pa, count_wrong_cons = trim_seqs(dict_aligned_genomes, min_per_taxa / 100.0, consensus / 100.0, min_perc_aa / 100.0)
            self.logger.info('Trimmed alignment from %d to %d AA (%d by minimum taxa percent, %d by consensus).' % (len(dict_aligned_genomes[dict_aligned_genomes.keys()[0]]),
                                                                                                                    len(trimmed_seqs[trimmed_seqs.keys()[0]]), count_wrong_pa, count_wrong_cons))
            self.logger.info('After trimming %d taxa have amino acids in <%.1f%% of columns.' % (
                len(pruned_seqs), min_perc_aa))

            if not os.path.exists(outdir):
                os.makedirs(outdir)
            seqfile_path = os.path.join(outdir, prefix + "_user_msa.fna")
            seqfile = open(seqfile_path, 'w')
            for key, value in trimmed_seqs.iteritems():
                seqfile.write(">{0}\n{1}\n".format(key, value))
            seqfile.close()

            if not notree:
                self.logger.info(
                    'Inferring tree for %d genomes under the WAG and GAMMA models.' % len(dict_aligned_genomes))

                output_tree = os.path.join(outdir, prefix + '_phylogeny.wag_gamma.tree')
                output_tree_log = os.path.join(outdir, 'fasttree.log')
                log_file = os.path.join(outdir, 'fasttree_output.txt')

                fasttree = FastTree(multithreaded=True)
                fasttree.run(
                    seqfile_path, 'prot', 'wag', output_tree, output_tree_log, log_file)

            self.logger.info('Done.')
            return True

        except IOError as e:
            self.logger.error(str(e))
            self.logger.error("GtdbTK has stopped.")

    def _selectGTDBGenomes(self, concatenated_file, filter_taxa):
        if filter_taxa is None:
            return read_fasta(concatenated_file)
        else:
            temp_dict = read_fasta(concatenated_file)
            list_required_fasta = filter_taxa.split(',')
            dict_taxonomy = {}
            with open(self.taxonomy_file, 'r') as taxfile:
                for line in taxfile:
                    tax_info = line.strip().split('\t')
                    genome_taxonomy = tax_info[1].split(';')
                    u = list(set(genome_taxonomy) & set(list_required_fasta))
                    if len(u) == 0:
                        if tax_info[0] in temp_dict:
                            del temp_dict[tax_info[0]]
            return temp_dict

    def list_genomes(self, batchfile):
        """Add genomes to database.

        Parameters
        ----------
        batchfile : str
            Name of file describing genomes to add.
        Returns
        -------
        bool
            True if genomes added without error.
        """
        # Add the genomes
        genomic_files = {}
        fh = open(batchfile, "rb")
        for line in fh:
            line = line.rstrip()
            if line == '':
                self.logger.warning(
                    "Encountered blank line in batch file. It has been ignored.")
                continue

            splitline = line.split("\t")
            if len(splitline) < 2:
                splitline += [None] * (2 - len(splitline))

            (fasta_path, name) = splitline

            if fasta_path is None or fasta_path == '':
                raise GenomeDatabaseError(
                    "Each line in the batch file must specify a path to the genome's fasta file.")

            if name is None or name == '':
                raise GenomeDatabaseError(
                    "Each line in the batch file must specify a name for the genome.")

            abs_fasta_path = os.path.abspath(fasta_path)
            real_fasta_path = os.path.realpath(fasta_path)

            if name in genomic_files and genome_files.get(name) == fasta_path:
                raise GenomeDatabaseError("{0} appears multiple times in the batchfile.".format(name))

            genomic_files[name] = fasta_path
        return genomic_files

    def _prepareGenomes(self, dict_user_genomes, output_dir):
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
            prodigal_output_dir = os.path.join(genome_output_dir, self.userAnnotationDir)
            if not os.path.exists(prodigal_output_dir):
                os.makedirs(prodigal_output_dir)

            # copy the genome file to the temporary folder
            fasta_target_file = os.path.join(genome_output_dir, user_genome + self.genomeFileSuffix)
            shutil.copy(dict_user_genomes.get(user_genome_raw), fasta_target_file)

            genomic_files[user_genome] = fasta_target_file
        return genomic_files

    def _alignMarkerSet(self, db_genome_ids, maindomain):
        manager = multiprocessing.Manager()
        out_q = manager.Queue()
        procs = []
        nprocs = self.threads
        for item in splitchunks(db_genome_ids, nprocs):
            p = multiprocessing.Process(
                target=self._hmmWorker,
                args=(item, out_q, maindomain))
            procs.append(p)
            p.start()

        # Collect all results into a single result dict. We know how many dicts
        # with results to expect.
        while out_q.empty():
            time.sleep(1)

        # Wait for all worker processes to finish
        results = {}

        for i in range(len(db_genome_ids)):
            id, sequence = out_q.get()
            results[id] = sequence

        return results

    def _calculatePresenceMarkerSets(self, db_genome_ids, outfile):
        """Generate a Tree based on the MSA.

        Parameters
        ----------
        user_directory : str
            Directory path where all fasta files are.
        tree_generate : boolean
            if True : a tree is created. if False, the domain of each genome is predicted
        outfile : str
            Path to the output file for domain prediction
        maindomain: str
            "bacteria" or "archaea". This flag informa which marker set and which MSA to use.
        taxa_filter: str
            Select the taxa of interest to generate the tree
        """

        manager = multiprocessing.Manager()
        out_q = manager.Queue()
        procs = []
        nprocs = self.threads
        for item in splitchunks(db_genome_ids, nprocs):
            p = multiprocessing.Process(
                target=self._hmmWorker,
                args=(item, out_q))
            procs.append(p)
            p.start()

        # Collect all results into a single result dict. We know how many dicts
        # with results to expect.
        while out_q.empty():
            time.sleep(1)

        # Wait for all worker processes to finish
        results = {"Bacteria": [], "Archaea": [], "Undefined": []}
        for i in range(len(db_genome_ids)):
            id, domain = out_q.get()
            if domain == "d__Bacteria":
                results.get("Bacteria").append(id)
            elif doamin == "d__Archaea":
                results.get("Archaea").append(id)
            else:
                results.get("Undefined").append(id)

        outf = open(outfile, "w")
        for domain in results:
            outf.write("##########\n")
            outf.write("{0}\n".format(domain))
            outf.write("##########\n")
            for genome in results[domain]:
                outf.write("{0}\n".format(genome))

        outf.close()
        return results

    def _hmmWorker(self, subdict_genomes, out_q, maindomain=None):
        '''
        The worker function, invoked in a process.

        :param subdict_genomes: sub dictionary of genomes
        :param out_q: manager.Queue()

        '''
        if maindomain:
            for db_genome_id, info in subdict_genomes.items():
                sequence = self._runHmmMultiAlign(db_genome_id, info.get("aa_gene_path"), maindomain)
                out_q.put((db_genome_id, sequence))
        else:
            for db_genome_id, info in subdict_genomes.items():
                domain = self._runHmmDetector(db_genome_id, info.get("aa_gene_path"))
                out_q.put((db_genome_id, domain))

        return True

    def _runHmmDetector(self, db_genome_id, path):
        '''
        Selects markers that are not present for a specific genome.

        :param db_genome_id: Selected genome
        :param path: Path to the genomic fasta file for the genome
        '''

        # gather information for all marker genes
        marker_dbs = {"PFAM": ConfigMetadata.PFAM_TOP_HIT_SUFFIX,
                      "TIGR": ConfigMetadata.TIGRFAM_TOP_HIT_SUFFIX}
        bacterial_set = Config.BACTERIAL_MARKERS
        archaeal_set = Config.ARCHAEAL_MARKERS

        marker_bac_dict_original, marker_arc_dict_original = [], []

        for db_marker in bacterial_set.keys():
            marker_bac_dict_original.extend([marker.replace(".HMM", "").replace(".hmm", "") for marker in bacterial_set[db_marker]])

        for db_marker in archaeal_set.keys():
            marker_arc_dict_original.extend([marker.replace(".HMM", "").replace(".hmm", "") for marker in archaeal_set[db_marker]])
        gene_bac_list, gene_arc_list = [], []

        for marker_db, marker_suffix in marker_dbs.iteritems():
            # get all gene sequences
            genome_path = str(path)
            tophit_path = genome_path.replace(ConfigMetadata.PROTEIN_FILE_SUFFIX, marker_suffix)

            # we store the tophit file line by line and store the
            # information in a dictionary
            with open(tophit_path) as tp:
                # first line is header line
                tp.readline()
                gene_dict = {}
                for line_tp in tp:
                    linelist = line_tp.split("\t")
                    genename = linelist[0]
                    sublist = linelist[1]
                    if ";" in sublist:
                        diff_markers = sublist.split(";")
                    else:
                        diff_markers = [sublist]

                    for each_gene in diff_markers:
                        sublist = each_gene.split(",")
                        markerid = sublist[0]
                        if markerid not in marker_bac_dict_original and markerid not in marker_arc_dict_original:
                            continue
                        if markerid in marker_bac_dict_original:
                            gene_bac_list.append(markerid)
                        if markerid in marker_arc_dict_original:
                            gene_arc_list.append(markerid)

        # gene_arc_list & gene_bac_list return the list of marker present in the genome
        # We can calculate the frequencies of markers for both Bacterial and Archaeal sets.
        arc_aa_per = (len(set(gene_arc_list)) * 100 / len(marker_arc_dict_original))
        bac_aa_per = (len(set(gene_bac_list)) * 100 / len(marker_bac_dict_original))
        # depending of the frequencies, a domain is assigned to the genome.
        if arc_aa_per < DefaultValues.DEFAULT_DOMAIN_THRESHOLD and bac_aa_per < DefaultValues.DEFAULT_DOMAIN_THRESHOLD:
            domain = "d__Undefined"
        elif bac_aa_per >= arc_aa_per:
            domain = "d__Bacteria"
        else:
            domain = "d__Archaea"

        return domain

    def _runHmmMultiAlign(self, db_genome_id, path, maindomain):
        '''
        Returns the concatenated marker sequence for a specific genome

        :param db_genome_id: Selected genome
        :param path: Path to the genomic fasta file for the genome
        :param maindomain: "bacteria" or "archaea" marker sets
        '''

        # gather information for all marker genes
        marker_dbs = {"PFAM": ConfigMetadata.PFAM_TOP_HIT_SUFFIX,
                      "TIGRFAM": ConfigMetadata.TIGRFAM_TOP_HIT_SUFFIX}

        marker_paths = {"PFAM": Config.PFAM_HMM_DIR,
                        "TIGRFAM": Config.TIGRFAM_HMM_DIR}

        bacterial_set = Config.BACTERIAL_MARKERS
        archaeal_set = Config.ARCHAEAL_MARKERS

        marker_dict_original = {}

        ordered_markers = []

        if maindomain == "bacteria":
            for db_marker in sorted(bacterial_set):
                marker_dict_original.update({marker.replace(".HMM", "").replace(".hmm", ""): os.path.join(marker_paths[db_marker], marker) for marker in bacterial_set[db_marker]})
        else:
            for db_marker in sorted(archaeal_set):
                marker_dict_original.update({marker.replace(".HMM", "").replace(".hmm", ""): os.path.join(marker_paths[db_marker], marker) for marker in archaeal_set[db_marker]})

        result_aligns = {}
        result_aligns[db_genome_id] = {}

        for marker_db, marker_suffix in marker_dbs.iteritems():
            # get all gene sequences
            genome_path = str(path)
            tophit_path = genome_path.replace(ConfigMetadata.PROTEIN_FILE_SUFFIX, marker_suffix)

            # we load the list of all the genes detected in the genome
            protein_file = tophit_path.replace(
                marker_suffix, ConfigMetadata.PROTEIN_FILE_SUFFIX)
            all_genes_dict = read_fasta(protein_file, False)
            # we store the tophit file line by line and store the
            # information in a dictionary
            with open(tophit_path) as tp:
                # first line is header line
                tp.readline()
                gene_dict = {}
                for line_tp in tp:
                    linelist = line_tp.split("\t")
                    genename = linelist[0]
                    sublist = linelist[1]
                    if ";" in sublist:
                        diff_markers = sublist.split(";")
                    else:
                        diff_markers = [sublist]

                    for each_gene in diff_markers:
                        sublist = each_gene.split(",")
                        markerid = sublist[0]
                        if markerid not in marker_dict_original.keys():
                            continue
                        evalue = sublist[1]
                        bitscore = sublist[2].strip()

                        if markerid in gene_dict:
                            oldbitscore = gene_dict.get(markerid).get("bitscore")
                            if oldbitscore < bitscore:
                                gene_dict[markerid] = {"marker_path": marker_dict_original.get(markerid),
                                                       "gene": genename,
                                                       "gene_seq": all_genes_dict.get(genename),
                                                       "bitscore": bitscore}
                        else:
                            gene_dict[markerid] = {"marker_path": marker_dict_original.get(markerid),
                                                   "gene": genename,
                                                   "gene_seq": all_genes_dict.get(genename),
                                                   "bitscore": bitscore}

            for mid, mpath in marker_dict_original.iteritems():
                if mid not in gene_dict and mid not in result_aligns.get(db_genome_id):
                    size = self._getHmmSize(mpath)
                    result_aligns.get(db_genome_id).update({mid: "-" * size})
                    #final_genome.append((db_genome_id, mid, "-" * size))

            result_aligns.get(db_genome_id).update(self._runHmmAlign(gene_dict, db_genome_id))
        # we concatenate the aligned markers together and associate them with the genome.
        for gid, markids in result_aligns.iteritems():
            seq = ""
            for markid in sorted(markids.keys()):
                seq = seq + markids.get(markid)

        return seq

    def _runHmmAlign(self, marker_dict, genome):
        '''
        Run hmmalign for a set of genes for a specific genome. This is run in a temp folder.

        :param marker_dict: list of markers that need to be aligned
        :param genome: specific genome id

        Returns
        --------------
        List of tuple to be inserted in aligned_markers table
        '''
        result_genomes_dict = {}
        hmmalign_dir = tempfile.mkdtemp()
        input_count = 0
        for markerid, marker_info in marker_dict.iteritems():
            hmmalign_gene_input = os.path.join(
                hmmalign_dir, "input_gene{0}.fa".format(input_count))
            input_count += 1
            out_fh = open(hmmalign_gene_input, 'wb')
            out_fh.write(">{0}\n".format(marker_info.get("gene")))
            out_fh.write("{0}".format(marker_info.get("gene_seq")))
            out_fh.close()
            proc = subprocess.Popen(["hmmalign", "--outformat", "Pfam", marker_info.get(
                "marker_path"), hmmalign_gene_input], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            proc.wait()

            for line in proc.stderr:
                print "TODO"
            result = self._getAlignedMarker(
                marker_info.get("gene"), proc.stdout)
            if len(result) < 1:
                return "TODO"
            result_genomes_dict[markerid] = result

            input_count += 1
        shutil.rmtree(hmmalign_dir)
        return result_genomes_dict

    def _getHmmSize(self, path):
        size = 0
        with open(path) as fp:
            for line in fp:
                if line.startswith("LENG  "):
                    size = line.split("  ")[1]
                    break
        return int(size)

    def _getAlignedMarker(self, hit_name, result_file):
        '''
        Parse the output of Hmmalign

        :param hit_name: gene name
        :param result_file: output file from Hmmalign
        '''
        hit_seq = None
        mask_seq = None

        for line in result_file:
            splitline = line.split(" ", 1)
            if splitline[0] == hit_name.split(" ", 1)[0]:
                rsplitline = line.rsplit(" ", 1)
                hit_seq = rsplitline[-1]
                continue
            if line[0:len("#=GC RF")] == "#=GC RF":
                rsplitline = line.rsplit(" ", 1)
                mask_seq = rsplitline[-1]

        if mask_seq is None:
            raise Exception("Unable to get mask from hmm align result file")

        if hit_seq is None:
            return None

        aligned_marker = ""
        for pos in xrange(0, len(mask_seq)):
            if mask_seq[pos] != 'x':
                continue
            aligned_marker += hit_seq[pos]
        return aligned_marker
