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
import multiprocessing
import time
import tempfile
import subprocess
import shutil

from biolib.external.execute import check_dependencies
from biolib.checksum import sha256
from biolib.seq_io import read_fasta

from ..tools import splitchunks, list_genomes_dir, merge_two_dicts


class HmmAligner(object):
    """Runs HMMalign over a set of genomes."""

    def __init__(self, 
                    threads,
                    pfam_top_hit_suffix,
                    tigrfam_top_hit_suffix,
                    protein_file_suffix,
                    pfam_hmm_dir,
                    tigrfam_hmm_dir,
                    bac120_markers,
                    ar122_markers,
                    rps23_markers):
        """Initialization."""
        
        check_dependencies(['hmmalign'])

        self.threads = threads
        self.pfam_top_hit_suffix = pfam_top_hit_suffix
        self.tigrfam_top_hit_suffix = tigrfam_top_hit_suffix
        self.protein_file_suffix = protein_file_suffix
        self.pfam_hmm_dir = pfam_hmm_dir
        self.tigrfam_hmm_dir = tigrfam_hmm_dir
        
        self.bac120_markers = bac120_markers
        self.ar122_markers = ar122_markers
        self.rps23_markers = rps23_markers
        
    def align_marker_set(self, db_genome_ids, marker_set_id):
        manager = multiprocessing.Manager()
        out_q = manager.Queue()
        procs = []
        nprocs = self.threads
        for item in splitchunks(db_genome_ids, nprocs):
            p = multiprocessing.Process(
                target=self._worker,
                args=(item, out_q, marker_set_id))
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
    
    def _worker(self, subdict_genomes, out_q, marker_set_id):
        '''
        The worker function, invoked in a process.
        :param subdict_genomes: sub dictionary of genomes
        :param out_q: manager.Queue()
        '''
        for db_genome_id, info in subdict_genomes.items():
            sequence = self._run_multi_align(db_genome_id, 
                                                info.get("aa_gene_path"), 
                                                marker_set_id)
            out_q.put((db_genome_id, sequence))
        return True
     
    def _run_multi_align(self, db_genome_id, path, marker_set_id):
        '''
        Returns the concatenated marker sequence for a specific genome
        :param db_genome_id: Selected genome
        :param path: Path to the genomic fasta file for the genome
        :param marker_set_id: Unique ID of marker set to use for alignment
        '''
        
        # gather information for all marker genes
        marker_paths = {"PFAM": os.path.join(self.pfam_hmm_dir,'individual_hmms'),
                        "TIGRFAM": os.path.join(os.path.dirname(self.tigrfam_hmm_dir),'individual_hmms')}
       
        marker_dict_original = {}
        if marker_set_id == "bac120":
            for db_marker in sorted(self.bac120_markers):
                marker_dict_original.update({marker.replace(".HMM", "").replace(".hmm", ""): 
                                                os.path.join(marker_paths[db_marker], marker) 
                                                for marker in self.bac120_markers[db_marker]})
        elif marker_set_id == "ar122":
            for db_marker in sorted(self.ar122_markers):
                marker_dict_original.update({marker.replace(".HMM", "").replace(".hmm", ""): 
                                                os.path.join(marker_paths[db_marker], marker) 
                                                for marker in self.ar122_markers[db_marker]})
        elif marker_set_id == "rps23":
            for db_marker in sorted(self.rps23_markers):
                marker_dict_original.update({marker.replace(".HMM", "").replace(".hmm", ""): 
                                                os.path.join(marker_paths[db_marker], marker) 
                                                for marker in self.rps23_markers[db_marker]})

        result_aligns = {}
        result_aligns[db_genome_id] = {}

        marker_dbs = {"PFAM": self.pfam_top_hit_suffix,
                      "TIGRFAM": self.tigrfam_top_hit_suffix}
        for marker_db, marker_suffix in marker_dbs.iteritems():
            # get all gene sequences
            genome_path = str(path)
            tophit_path = genome_path.replace(self.protein_file_suffix, marker_suffix)

            # we load the list of all the genes detected in the genome
            protein_file = tophit_path.replace(
                marker_suffix, self.protein_file_suffix)
            all_genes_dict = read_fasta(protein_file, False)
            
            # Prodigal adds an asterisks at the end of each called genes,
            # These asterisks sometimes appear in the MSA, which can be an issue for some softwares downstream
            for seq_id, seq in all_genes_dict.iteritems():
                if seq[-1] == '*':
                    all_genes_dict[seq_id] = seq[:-1]
            
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
                    size = self._get_hmm_size(mpath)
                    result_aligns.get(db_genome_id).update({mid: "-" * size})
                    #final_genome.append((db_genome_id, mid, "-" * size))

            result_aligns.get(db_genome_id).update(self._run_align(gene_dict, db_genome_id))
        
        # we concatenate the aligned markers together and associate them with the genome.
        for gid, markids in result_aligns.iteritems():
            seq = ""
            for markid in sorted(markids.keys()):
                seq = seq + markids.get(markid)

        return seq
    
    def _run_align(self, marker_dict, genome):
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
                print line
                
            result = self._get_aligned_marker(marker_info.get("gene"), proc.stdout)                
            if len(result) < 1:
                return "TODO"
                
            result_genomes_dict[markerid] = result

            input_count += 1
        shutil.rmtree(hmmalign_dir)
        return result_genomes_dict

    def _get_hmm_size(self, path):
        size = 0
        with open(path) as fp:
            for line in fp:
                if line.startswith("LENG  "):
                    size = line.split("  ")[1]
                    break
        return int(size)
    
    def _get_aligned_marker(self, hit_name, result_file):
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
    