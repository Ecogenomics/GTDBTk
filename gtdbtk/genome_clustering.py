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
from statistics import mean

from gtdbtk.biolib_lite.seq_io import read_fasta
from gtdbtk.config.common import CONFIG
from gtdbtk.exceptions import GTDBTkExit
from gtdbtk.external.fastani import FastANI


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
        self.contig_count = len(seqs)
        self.calculateGenomeScore()

    def calculateGenomeScore(self):
        """ Calculate the genome score based on the completeness and contamination values ,number of contigs and ambiguous bases."""
        formula = self.completeness - 5 * self.contamination
        formula = formula - 5 * (self.contig_count/100)
        formula = formula - 5 * (self.ambiguous_bases/10000)
        self.score = formula


class GenomeClustering(object):
    """ Store the information from the CheckM summary file.
     The template should be able to handle both version 1 and 2 of CheckM but also generic completeness and contamination files"""

    checkm2_headers =['Name','Completeness_General','Contamination','Contig_N50']
    checkm1_headers =['Bin Id','Completeness','Contamination']
    generic_headers = ['Name','Completeness','Contamination']

    def __init__(self,path:str,genome_path:dict,cpus:int = 1):
        self.logger = logging.getLogger('timestamp')
        self.genome_path = genome_path
        self.de_novo_file_path: str = path
        self.de_novo_file_format=0# 0: generic, 1: checkm1, 2: checkm2
        self.qc_genomes = None
        self.cpus=cpus

    def parse(self,genomes_to_process):
        genomes = dict()
        name_index, completeness_index, contamination_index,sep = self.__detect_format_checkm()

        with open(self.de_novo_file_path, 'r') as fh:
            fh.readline()
            for line in fh:
                line = line.strip().split(sep)
                filepath = self.genome_path[line[name_index]]
                if line[name_index] not in self.genome_path:
                    raise GTDBTkExit(f"Genome {line[name_index]} not found in the genome folder/batch file.")
                if genomes_to_process and line[name_index] in genomes_to_process:
                    qg = QCedGenome(line[name_index],filepath, float(line[completeness_index]), float(line[contamination_index]))
                    qg.calculateSeqStats()
                    genomes[qg.name] = qg.score
        self.qc_genomes=genomes
        return genomes

    def species_cluster(self):
        fastani = FastANI(self.cpus, force_single=True)
        qc_genomes = self.qc_genomes
        #we sort the genomes by score
        sorted_genomes = sorted(qc_genomes.items(), key=lambda x: x[1], reverse=True)
        #we pick the first one as rep
        rep_id = sorted_genomes[0][0]
        #we compare the rep to all the other genomes
        other_genomes = list(zip(*sorted_genomes[1:]))[0]
        # we run fastANI between the rep and all the other genomes
        # we keep the ones that are above 95% identity
        temp_dict = {rep_id:other_genomes}
        fastani_results = fastani.run(temp_dict, self.genome_path)
        for gid in fastani_results.keys():
            thresh_results = [(ref_gid, hit) for (ref_gid, hit) in fastani_results[gid].items() if
                              hit['af'] >= CONFIG.AF_THRESHOLD and hit['ani'] >= CONFIG.FASTANI_SPECIES_THRESHOLD]
            closest = sorted(thresh_results, key=lambda x: (-x[1]['ani'], -x[1]['af']))
            if len(closest) > 0:
                a=1


    def __detect_format_checkm(self):
        with open(self.de_novo_file_path) as f:
            header = f.readline()
        sep = self.__detect_delimiters(header)
        header_infos = header.strip().split(sep)
        if all(item in header_infos for item in self.checkm2_headers):
            self.format = 2
            name_index = header_infos.index('Name')
            completeness_index = header_infos.index('Completeness_General')
            contamination_index = header_infos.index('Contamination')
            return name_index, completeness_index, contamination_index ,sep
        elif all(item in header_infos for item in self.checkm1_headers):
            self.format = 1
            name_index = header_infos.index('Bin Id')
            completeness_index = header_infos.index('Completeness')
            contamination_index = header_infos.index('Contamination')
            return name_index, completeness_index, contamination_index,sep
        elif all(item in header_infos for item in self.generic_headers):
            self.format = 0
            name_index = header_infos.index('Name')
            completeness_index = header_infos.index('Completeness')
            contamination_index = header_infos.index('Contamination')
            return name_index, completeness_index, contamination_index, sep
        else:
            raise Exception('Unknown format of Completeness/Contamination file. the format should be either CheckM1, CheckM2 '
                            'or contains the headers Name, Completeness and Contamination.')

    def __detect_delimiters(self,header):
        delimiters = ['\t', ',']
        for delimiter in delimiters:
            if delimiter in header:
                return delimiter
        return None

