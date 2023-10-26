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
import shutil
from collections import defaultdict
from typing import List, Tuple

from gtdbtk.biolib_lite.common import make_sure_path_exists
from gtdbtk.biolib_lite.seq_io import read_fasta
from gtdbtk.config.common import CONFIG
from gtdbtk.config.output import DIR_ANI_REP_INT_DEREP, DIR_ANI_REP_INT_MASH, DIR_ANISCREEN
from gtdbtk.exceptions import GTDBTkExit
from gtdbtk.external.fastani import FastANI
from gtdbtk.external.mash import Mash
from gtdbtk.files.classify_summary import ClassifySummaryFileRow
from gtdbtk.tools import tqdm_log


class QCedGenome:
    """ Class for individual genome quality metrics"""

    CONTIG_BREAK = 'NNNNNNNNNN'

    def __init__(self,name,path,completeness,contamination):
        self.name = name
        self.path = path
        self.completeness = completeness
        self.contamination = contamination
        self.contig_count = None
        self.ambiguous_bases = None
        self.score = None


    def calculateSeqStats(self):
        """ Based on CheckM code. Calculate scaffold length statistics (min length, max length, total length, N50, # contigs)."""

        # read scaffolds
        seqs = read_fasta(self.path)
        numAmbiguousBases = 0
        for seqId, seq in seqs.items():

            numAmbiguousBases += seq.count('N') + seq.count('n')
        self.ambiguous_bases = numAmbiguousBases
        self.contig_count = len(self.identify_contigs(seqs))
        self.calculateGenomeScore()

    def identify_contigs(self,seqs, contig_break=CONTIG_BREAK):
        """Break scaffolds into contigs.

        Parameters
        ----------
        seqs : dict[seq_id] -> seq
            Sequences indexed by sequence ids.
        contig_break : str
            Motif used to split scaffolds into contigs.

        Returns
        -------
        dict : dict[seq_id] -> seq
            Contigs indexed by sequence ids.
        """

        contigs = {}
        for seq_id, seq in seqs.items():
            seq = seq.upper()
            contig_count = 0
            for contig in seq.split(contig_break):
                contig = contig.strip('N')
                if contig:
                    contigs[seq_id + '_c' + str(contig_count)] = contig
                    contig_count += 1

        return contigs

    def calculateGenomeScore(self):
        """ Calculate the genome score based on the completeness and contamination values ,number of contigs and ambiguous bases."""
        formula = self.completeness - 5 * self.contamination
        formula = formula - 5 * (self.contig_count/100)
        formula = formula - 5 * (self.ambiguous_bases/10000)
        self.score = formula


class GenomeClustering(object):
    """ Store the information from the CheckM summary file.
     The template should be able to handle both version 1 and 2 of CheckM but also generic completeness and contamination files"""

    checkm2_headers_allmodels = ['Name','Completeness_Model_Used']
    checkm2_headers_specific='Completeness_Specific'
    checkm2_headers_general= 'Completeness_General'
    checkm1_headers =['Bin Id','Completeness','Contamination']
    generic_headers = ['Name','Completeness','Contamination']

    def __init__(self,path:str,genome_path:dict,user_genomes_radii:dict,output_dir,prefix,max_d, mash_k, mash_v, mash_s,cpus:int = 1):
        self.logger = logging.getLogger('timestamp')
        self.genome_path = genome_path
        self.de_novo_file_path: str = path
        self.de_novo_file_format=0# 0: generic, 1: checkm1, 2: checkm2
        self.user_genomes_radii = user_genomes_radii
        self.qc_genomes = None
        self.cpus=cpus
        self.intermediates_directory = output_dir
        self.prefix = prefix
        self.max_d =max_d
        self.mash_k = mash_k
        self.mash_v = mash_v
        self.mash_s = mash_s


    def parse(self,genomes_to_process):
        genomes = dict()
        name_index, completeness_index, contamination_index,sep,both_model = self.__detect_format_checkm()
        with open(self.de_novo_file_path, 'r') as fh:
            fh.readline()
            for line in fh:
                line = line.strip().split(sep)
                filepath = self.genome_path[line[name_index]]
                if line[name_index] not in self.genome_path:
                    raise GTDBTkExit(f"Genome {line[name_index]} not found in the genome folder/batch file.")
                if genomes_to_process and line[name_index] in genomes_to_process:
                    if both_model:
                        # we need to find the right completeness value,completeness_index needs to be parsed
                        model_info_idx = both_model.get('model_used_idx')
                        if '(General Model)' in line[model_info_idx]:
                            completeness_index = both_model.get('general_idx')
                        elif '(Specific Model)' in line[model_info_idx]:
                            completeness_index = both_model.get('specific_idx')
                    qg = QCedGenome(line[name_index],filepath, float(line[completeness_index]), float(line[contamination_index]))
                    qg.calculateSeqStats()
                    genomes[qg.name] = qg.score

        self.qc_genomes=genomes
        return genomes

    def species_cluster(self,all_genomes_to_process,marker_set_id):
        fastani = FastANI(self.cpus, force_single=True)
        qc_genomes = self.qc_genomes
        #we sort the genomes by score
        sorted_genomes = sorted(qc_genomes.items(), key=lambda x: x[1], reverse=True)
        len_sorted_genomes = len(sorted_genomes)
        sp_cluster= defaultdict()
        increment = 1

        dir_mash = os.path.join(self.intermediates_directory, DIR_ANI_REP_INT_DEREP,'cluster',marker_set_id)
        mash = Mash(self.cpus, dir_mash, self.prefix)

        ref_genomes = {gid: self.genome_path[gid] for gid in all_genomes_to_process}
        mash_results = mash.run(ref_genomes, ref_genomes, self.max_d, self.mash_k, self.mash_v, self.mash_s, mash_max_dist=100, mash_db=dir_mash,all_vs_all = True)

        with tqdm_log(unit='genome', total=len_sorted_genomes) as p_bar:
            while len(sorted_genomes) > 1:
                #we pick the first one as rep
                rep_id = sorted_genomes[0][0]
                sp_cluster[rep_id] = {'increment':increment,'members':[],'radius':self.user_genomes_radii[rep_id]}
                #we compare the rep to all the other genomes
                other_genomes = list(zip(*sorted_genomes[1:]))[0]
                # we run fastANI between the rep and all the other genomes
                # we keep the ones that are above 95% identity
                # we can remove all genomes that are not in the pre filter from mash
                other_genomes = [gid for gid in other_genomes if gid in mash_results[rep_id]]
                if other_genomes:
                    temp_dict = {rep_id:other_genomes}
                    fastani_results = fastani.run(temp_dict, self.genome_path, verbose=False)
                else:
                    fastani_results = dict()
                processed_genomes = (rep_id,)

                if rep_id in fastani_results.keys():
                    thresh_results = [(ref_gid, hit) for (ref_gid, hit) in fastani_results[rep_id].items() if
                                      hit['af'] >= CONFIG.AF_THRESHOLD
                                      and hit['ani'] >= (max(self.user_genomes_radii[rep_id],self.user_genomes_radii[ref_gid]))]
                    closest = sorted(thresh_results, key=lambda x: (-x[1]['ani'], -x[1]['af']))
                    sp_cluster[rep_id]['members']=closest
                    if len(closest)>0:
                        processed_genomes = (rep_id,) + list(zip(*closest))[0]
                sorted_genomes = [(i, v) for i, v in sorted_genomes if i not in processed_genomes]
                increment += 1
            #we add the last rep
                p_bar.update(len(processed_genomes))
            if len(sorted_genomes) == 1:
                rep_id = sorted_genomes[0][0]
                sp_cluster[rep_id] = {'increment':increment,'members':[]}
                p_bar.update(1)

        return sp_cluster

    def compare_gtdb_reps(self,qc_genomes_cluster: dict ,marker_set_id,  list_summary_rows: List[ClassifySummaryFileRow]) -> Tuple[dict, dict]:
        fastani = FastANI(self.cpus, force_single=True)
        value_ani = defaultdict()
        changing_reps_genomes = []
        # we run mash before fastANI
        dir_mash = os.path.join(self.intermediates_directory, DIR_ANI_REP_INT_DEREP, 'compare', marker_set_id)
        make_sure_path_exists(dir_mash)
        mash = Mash(self.cpus, dir_mash, self.prefix)

        ref_genomes = {gid: self.genome_path[gid] for gid,clusterg in qc_genomes_cluster.items()}
        genomes_to_process = {}
        for row in list_summary_rows:
            genomes_to_process[row.gid] = self.genome_path[row.gid]
            value_ani[row.gid] = {'rep': row.fastani_ref, 'ani': row.fastani_ani, 'af': row.fastani_af}
        for repg,clustered_gs in qc_genomes_cluster.items():
            for clustered_g in clustered_gs['members']:
                genomes_to_process[clustered_g[0]] = self.genome_path[clustered_g[0]]
                value_ani[clustered_g[0]] = {'rep': repg, 'ani': clustered_g[1].get('ani'),
                                             'af': clustered_g[1].get('af')}

        mash_results = mash.run(genomes_to_process, ref_genomes, self.max_d, self.mash_k, self.mash_v, self.mash_s, mash_max_dist=100, mash_db=None)
        d_compare = defaultdict(set)
        for qry_gid, ref_hits in mash_results.items():
            d_compare[qry_gid] = d_compare[qry_gid].union(set(ref_hits.keys()))
        for row in list_summary_rows:
            value_ani[row.gid] = {'rep':row.fastani_ref,'ani':row.fastani_ani,'af':row.fastani_af}
        for repg,clustered_gs in qc_genomes_cluster.items():
            for clustered_g in clustered_gs['members']:
                value_ani[clustered_g[0]] = {'rep':repg,'ani':clustered_g[1].get('ani'),'af':clustered_g[1].get('af')}
        fastani_results = fastani.run(d_compare, self.genome_path)
        for gid in fastani_results.keys():
            thresh_results = [(ref_gid, hit) for (ref_gid, hit) in fastani_results[gid].items() if
                              hit['af'] >= CONFIG.AF_THRESHOLD and hit['ani'] >= CONFIG.FASTANI_SPECIES_THRESHOLD]
            closest = sorted(thresh_results, key=lambda x: (-x[1]['ani'], -x[1]['af']))
            if closest:
                closest_user_rep = closest[0]
                if closest_user_rep[0] != value_ani[gid]['rep']:
                    if closest_user_rep[1]['ani'] > value_ani[gid]['ani'] and closest_user_rep[1]['af'] >= value_ani[gid]['af']:
                        qc_genomes_cluster[closest_user_rep[0]]['members'].append(gid)
                        changing_reps_genomes.append(gid)
                    elif closest_user_rep[1]['ani'] >= value_ani[gid]['ani'] and closest_user_rep[1]['af'] > value_ani[gid]['af']:
                        qc_genomes_cluster[closest_user_rep[0]]['members'].append(gid)
                        changing_reps_genomes.append(gid)


        return qc_genomes_cluster, changing_reps_genomes


    def __detect_format_checkm(self):
        with open(self.de_novo_file_path) as f:
            header = f.readline()
        sep = self.__detect_delimiters(header)
        both_models = {}
        header_infos = header.strip().split(sep)
        if all(item in header_infos for item in self.checkm2_headers_allmodels) and self.checkm2_headers_general in header_infos and self.checkm2_headers_specific in header_infos:
            self.logger.info('CheckM 2 file format detected')
            name_index = header_infos.index('Name')
            completeness_index = None
            contamination_index = header_infos.index('Contamination')
            both_models = {'model_used_idx': header_infos.index('Completeness_Model_Used'),
                           'general_idx': header_infos.index(self.checkm2_headers_general),
                           'specific_idx': header_infos.index(self.checkm2_headers_specific)}
            return name_index, completeness_index, contamination_index, sep, both_models
        elif self.checkm2_headers_general in header_infos and self.checkm2_headers_specific not in header_infos:
            self.de_novo_file_format = 2
            self.logger.info('CheckM 2 file format detected')
            name_index = header_infos.index('Name')
            completeness_index = header_infos.index(self.checkm2_headers_general)
            contamination_index = header_infos.index('Contamination')
            return name_index, completeness_index, contamination_index ,sep, both_models
        elif self.checkm2_headers_general not in header_infos and self.checkm2_headers_specific in header_infos:
            self.de_novo_file_format = 2
            self.logger.info('CheckM 2 file format detected')
            name_index = header_infos.index('Name')
            completeness_index = header_infos.index(self.checkm2_headers_specific)
            contamination_index = header_infos.index('Contamination')
            return name_index, completeness_index, contamination_index ,sep, both_models
        elif all(item in header_infos for item in self.checkm1_headers):
            self.logger.info('CheckM 1 file format detected')
            self.de_novo_file_format = 1
            name_index = header_infos.index('Bin Id')
            completeness_index = header_infos.index('Completeness')
            contamination_index = header_infos.index('Contamination')
            return name_index, completeness_index, contamination_index,sep, both_models
        elif all(item in header_infos for item in self.generic_headers):
            self.logger.info('Generic file format detected')
            self.de_novo_file_format = 0
            name_index = header_infos.index('Name')
            completeness_index = header_infos.index('Completeness')
            contamination_index = header_infos.index('Contamination')
            return name_index, completeness_index, contamination_index, sep, both_models
        else:
            raise Exception('Unknown format of Completeness/Contamination file. the format should be either CheckM1, CheckM2 '
                            'or contains the headers Name, Completeness and Contamination.')

    def __detect_delimiters(self,header):
        delimiters = ['\t', ',']
        for delimiter in delimiters:
            if delimiter in header:
                return delimiter
        return None

