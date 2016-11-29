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

from biolib.seq_io import read_fasta

from Tools import fastaPathGenerator, splitchunks, list_genomes_dir


class TreeManager(object):
    """Manages genomes, concatenated alignment, and metadata for tree inference and visualization."""

    def __init__(self):
        self.init = True
        self.logger = logging.getLogger()
        self.logger.setLevel(logging.DEBUG)

        ch = logging.StreamHandler(sys.stdout)
        ch.setLevel(logging.DEBUG)
        formatter = logging.Formatter('%(asctime)s - <%(levelname)s>: %(message)s')
        ch.setFormatter(formatter)
        self.logger.addHandler(ch)

    def MakeTreeData(self, user_directory, maindomain, taxa_filter=None):
        # create tree data
        try:
            dictionary_user_genomes = list_genomes_dir(user_directory)
            self.logger.info('Creating tree data for {0} genomes.'.format(len(dictionary_user_genomes)))

            genomic_files = self._addGenomeBatch(dict_user_genomes, self.tmp_output_dir)

            self.logger.info("Running Prodigal to identify genes.")
            prodigal = Prodigal(self.threads)
            genome_dictionary = prodigal.run(genomic_files)

            # annotated genes against TIGRfam and Pfam databases
            self.logger.info("Identifying TIGRfam protein families.")
            gene_files = [genome_dictionary[db_genome_id]['aa_gene_path']
                          for db_genome_id in genome_dictionary.keys()]
            tigr_search = TigrfamSearch(self.threads)
            tigr_search.run(gene_files)

            self.logger.info("Identifying Pfam protein families.")
            pfam_search = PfamSearch(self.threads)
            pfam_search.run(gene_files)

            dict_aligned_genomes = self._alignMarkerSet(genome_dictionary, maindomain)

            self.logger.info('Done.')
        except IOError as e:
            self.logger.error(str(e))
            self.logger.error("GtdbTK has stopped.")

        return True

    def _alignMarkerSet(self, genome_dictionary, maindomain):
        manager = multiprocessing.Manager()
        out_q = manager.Queue()
        procs = []
        nprocs = self.threads
        print db_genome_ids
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
        while not out_q.empty():
            id, domain = out_q.get()
            if domain == "d__Bacteria":
                results.get("Bacteria").append(id)
            elif doamin == "d__Archaea":
                results.get("Archaea").append(id)
            else:
                results.get("Undefined").append(id)

        return results

    def _hmmWorker(self, subdict_genomes, out_q, maindomain):
        '''
        The worker function, invoked in a process.

        :param subdict_genomes: sub dictionary of genomes
        :param out_q: manager.Queue()

        '''

        for db_genome_id, info in subdict_genomes.items():
            domain = self._runHmmMultiAlign(db_genome_id, info.get("aa_gene_path"), maindomain)
        out_q.put((db_genome_id, domain))
        return True

    def _runHmmMultiAlign(self, db_genome_id, path, maindomain):
        '''
        Selects markers that are not aligned for a specific genome.

        :param db_genome_id: Selected genome
        :param path: Path to the genomic fasta file for the genome
        :param marker_ids: list of marker ids for the selected sets
        '''

        # gather information for all marker genes
        marker_dbs = {"PFAM": ConfigMetadata.PFAM_TOP_HIT_SUFFIX,
                      "TIGR": ConfigMetadata.TIGRFAM_TOP_HIT_SUFFIX}

        marker_paths = {"PFAM": Config.PFAM_HMM_DIR,
                        "TIGR": Config.TIGRFAM_HMM_DIR}

        bacterial_set = Config.BACTERIAL_MARKERS
        archaeal_set = Config.ARCHAEAL_MARKERS

        marker_dict_original = []

        final_genome = []

        if maindomain == "bacteria":
            for db_marker in bacterial_set.keys().sort():
                marker_dict_original = {marker.replace(".HMM", "").replace(".hmm", ""): os.path.join(marker_paths[db_marker], marker) for marker in bacterial_set[db_marker]}
        else:
            for db_marker in archaeal_set.keys().sort():
                marker_dict_original = {marker.replace(".HMM", "").replace(".hmm", ""): os.path.join(marker_paths[db_marker], marker) for marker in archaeal_set[db_marker]}

        gene_bac_list = []
        gene_arc_list = []
        for marker_db, marker_suffix in marker_dbs.iteritems():
            # get all gene sequences
            genome_path = str(path)
            tophit_path = genome_path.replace(ConfigMetadata.PROTEIN_FILE_SUFFIX, marker_suffix)

            # we load the list of all the genes detected in the genome
            protein_file = tophit_path.replace(
                marker_suffix, self.protein_file_suffix)
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
                        if markerid not in marker_dict_original:
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

            for mid, path in marker_dict_original.iteritems():
                if mid not in gene_dict:
                    size = self._getHmmSize(path)
                    final_genome.append(mid, "-" * size)

            result_aligns = self._runHmmAlign(gene_dict, db_genome_id)
            for result_align in result_aligns:
                final_genome.append((result_align[0], result_align[2]))

        print final_genome
        return domain

    def _getHmmSize(self, path):
        size = 0
        with open(path) as fp:
            for line in fp:
                if line.startswith("LENG  "):
                    size = line.split("  ")[1]
                    break
        return size
