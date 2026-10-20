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
from typing import Dict, List, Tuple

from gtdbtk.biolib_lite.common import canonical_gid
from gtdbtk.biolib_lite.seq_io import read_fasta
from gtdbtk.config.common import CONFIG
from gtdbtk.exceptions import GTDBTkExit, GenomeNotInQCFile
from gtdbtk.external.skani import SkANI
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
        self.calculate_seq_stats()

    def calculate_seq_stats(self):
        """ Based on CheckM code. Calculate the number of contigs and
        ambiguous bases, then the genome score."""

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

    def can_be_representative(self):
        """ A genome can become a de novo representative only if it passes the quality gate."""
        return self.completeness > CONFIG.MIN_REP_COMPLETENESS and self.contamination < CONFIG.MAX_REP_CONTAMINATION


class GenomeClustering(object):
    """ Cluster user genomes into de novo species.

    The rules are the same as the species assignment in GTDB-Tk 2.8:
        1. a genome goes to its closest representative (AF >= af_threshold),
           and only the radius of that representative is checked;
        2. the radius of a new representative is max(95, ANI to its closest
           representative, GTDB or de novo), and is fixed once set;
        3. GTDB radii come from the release and never change;
        4. representatives are picked by decreasing quality score, and only
           genomes passing the quality gate can become representatives.

    The steps are:
        species_cluster    pick the de novo representatives among the genomes
                           not assigned to a GTDB species;
        recluster_genomes  final pass: every user genome goes to its closest
                           representative (GTDB or de novo) if it is within
                           its radius;
        cluster            runs both steps, and sends the genomes left without
                           a species back to species_cluster until none of them
                           can become a representative.

    ANI results are stored as ani_results[genome][other_genome] = {'ani': .., 'af': ..},
    with user genome pairs in both directions.
     """

    checkm2_headers_allmodels = ['Name','Completeness_Model_Used']
    checkm2_headers_specific='Completeness_Specific'
    checkm2_headers_general= 'Completeness_General'
    checkm1_headers =['Bin Id','Completeness','Contamination']
    generic_headers = ['Name','Completeness','Contamination']

    def __init__(self,path:str,genome_path:dict,gtdb_radii,output_dir,prefix,cpus:int = 1,
                 af_threshold:float = CONFIG.AF_THRESHOLD):
        self.logger = logging.getLogger('timestamp')
        self.genome_path = genome_path
        self.de_novo_file_path: str = path
        self.de_novo_file_format=0# 0: generic, 1: checkm1, 2: checkm2
        self.gtdb_radii = gtdb_radii
        self.qc_genomes = None
        self.cpus=cpus
        self.intermediates_directory = output_dir
        self.prefix = prefix
        self.af_threshold = af_threshold
        self.min_radius = CONFIG.SKANI_SPECIES_THRESHOLD
        # user genome ids, set by parse(); every other genome id is a GTDB representative
        self.user_genomes = set()

    def parse(self,genomes_to_process=None):
        """ Read the QC file for the genomes to process (all genomes if None).
        Rows for other genomes are ignored."""
        if genomes_to_process is None:
            genomes_to_process = set(self.genome_path.keys())
        genomes_to_process = set(genomes_to_process)
        genomes = dict()
        name_index, completeness_index, contamination_index,sep,both_model = self.__detect_format_checkm()
        with open(self.de_novo_file_path, 'r') as fh:
            fh.readline()
            for line in fh:
                if not line.strip():
                    continue
                line = line.rstrip('\n').split(sep)
                genome_id = line[name_index].strip()
                if genome_id not in genomes_to_process:
                    continue
                if genome_id not in self.genome_path:
                    raise GTDBTkExit(f"Genome {genome_id} not found in the genome folder/batch file.")
                if both_model:
                    # we need to find the right completeness value,completeness_index needs to be parsed
                    model_info_idx = both_model.get('model_used_idx')
                    if '(General Model)' in line[model_info_idx]:
                        completeness_index = both_model.get('general_idx')
                    elif '(Specific Model)' in line[model_info_idx]:
                        completeness_index = both_model.get('specific_idx')
                    else:
                        raise GTDBTkExit(f'Unknown CheckM2 model for {genome_id}: {line[model_info_idx]}')
                qg = QCedGenome(genome_id,self.genome_path[genome_id], float(line[completeness_index]), float(line[contamination_index]))
                genomes[qg.name] = qg

        missing_genomes = sorted(genomes_to_process - set(genomes.keys()))
        if missing_genomes:
            raise GenomeNotInQCFile(f'{len(missing_genomes)} genome(s) are not in the de_novo_species file, '
                                    f'e.g. {missing_genomes[0]}')

        self.qc_genomes=genomes
        self.user_genomes = set(genomes.keys())
        return {gid: qg.score for gid, qg in genomes.items()}

    def compute_ani(self, genome_ids):
        """ Run skani triangle (all-vs-all) between user genomes.
        Returns ani_results with both directions."""
        genomes = {gid: self.genome_path[gid] for gid in sorted(genome_ids)}
        ani_results = dict()
        if len(genomes) < 2:
            return ani_results
        self.logger.info(f'Running skani all-vs-all on {len(genomes):,} genomes.')
        skani = SkANI(self.cpus, force_single=True)
        skani_results = skani.run_triangle(genomes, self.prefix)
        for qry, hits in skani_results.items():
            for ref, hit in hits.items():
                self.add_hit(ani_results, qry, ref, hit)
        return ani_results

    @staticmethod
    def add_hit(ani_results, gid_a, gid_b, hit):
        """ Store a hit in both directions. If the pair is already there, the
        higher (ANI, AF) is kept so the result does not depend on the order."""
        if gid_a == gid_b:
            return
        current = ani_results.get(gid_a, {}).get(gid_b)
        if current is not None and (current['ani'], current['af']) >= (hit['ani'], hit['af']):
            return
        ani_results.setdefault(gid_a, {})[gid_b] = hit
        ani_results.setdefault(gid_b, {})[gid_a] = hit

    def get_radius(self, rep_id, sp_cluster):
        """ Radius of a de novo representative, or of a GTDB representative (from the release)."""
        if rep_id in sp_cluster:
            return sp_cluster[rep_id]['radius']
        return self.gtdb_radii.get_rep_ani(canonical_gid(rep_id))

    def closest_representative(self, gid, ani_results, sp_cluster):
        """ Closest representative (GTDB or de novo) passing the AF threshold:
        highest ANI, then AF, then genome id. Returns (rep_id, hit) or None."""
        # candidates: GTDB representatives and de novo representatives, never plain user genomes
        candidates = [(ref_id, hit) for ref_id, hit in ani_results.get(gid, {}).items()
                      if ref_id != gid
                      and (ref_id in sp_cluster or ref_id not in self.user_genomes)
                      and hit['af'] >= self.af_threshold]
        if not candidates:
            return None
        return sorted(candidates, key=lambda x: (-x[1]['ani'], -x[1]['af'], x[0]))[0]

    def is_within_radius(self, closest, sp_cluster):
        return closest is not None and closest[1]['ani'] >= self.get_radius(closest[0], sp_cluster)

    def species_cluster(self, genomes_to_cluster, ani_results, sp_cluster=None):
        """ Pick new representatives among genomes_to_cluster, best score first.
        A genome becomes a representative if its closest representative does not contain it."""
        if sp_cluster is None:
            sp_cluster = dict()
        increment = len(sp_cluster) + 1
        # best score first; the genome id breaks ties so the QC file order does not matter
        sorted_genomes = sorted(genomes_to_cluster, key=lambda gid: (-self.qc_genomes[gid].score, gid))

        with tqdm_log(unit='genome', total=len(sorted_genomes)) as p_bar:
            for gid in sorted_genomes:
                p_bar.update(1)
                if gid in sp_cluster or not self.qc_genomes[gid].can_be_representative():
                    continue
                closest = self.closest_representative(gid, ani_results, sp_cluster)
                if self.is_within_radius(closest, sp_cluster):
                    continue
                radius = self.min_radius if closest is None else max(self.min_radius, closest[1]['ani'])
                sp_cluster[gid] = {'increment': increment, 'members': [], 'radius': radius}
                increment += 1
        return sp_cluster

    def recluster_genomes(self, sp_cluster, gtdb_assignment, ani_results):
        """ Final pass: every user genome that is not a representative goes to its
        closest representative (GTDB or de novo) if it is within its radius.

        Returns
        -------
        sp_cluster : the de novo clusters, with their members rebuilt from scratch.
        gtdb_members : dict[gid] = (gtdb_rep, hit) for genomes in a GTDB species.
        changing_reps_genomes : dict[gid] = (old gtdb_rep, new rep or None) for genomes
            assigned to a GTDB species by the ANI screen whose species changed.
        unassigned_genomes : genomes with no species.
        """
        for rep_id in sp_cluster:
            sp_cluster[rep_id]['members'] = []
        gtdb_members = dict()
        changing_reps_genomes = dict()
        unassigned_genomes = set()

        for gid in sorted(self.user_genomes):
            # representatives are never re-assigned (a later rep R can have radius == ANI(R, this rep))
            if gid in sp_cluster:
                if gid in gtdb_assignment:
                    changing_reps_genomes[gid] = (gtdb_assignment[gid][0], gid)
                continue
            closest = self.closest_representative(gid, ani_results, sp_cluster)
            if self.is_within_radius(closest, sp_cluster):
                if closest[0] in sp_cluster:
                    sp_cluster[closest[0]]['members'].append((gid, closest[1]))
                else:
                    gtdb_members[gid] = closest
                new_rep = closest[0]
            else:
                unassigned_genomes.add(gid)
                new_rep = None
            if gid in gtdb_assignment and gtdb_assignment[gid][0] != new_rep:
                changing_reps_genomes[gid] = (gtdb_assignment[gid][0], new_rep)

        return sp_cluster, gtdb_members, changing_reps_genomes, unassigned_genomes

    def cluster(self, gtdb_assignment, ani_results):
        """ Run species_cluster and recluster_genomes until no genome left
        without a species can become a representative.

        Parameters
        ----------
        gtdb_assignment : dict[gid] = (gtdb_rep, hit)
            Genomes assigned to a GTDB species by the ANI screen.
        ani_results : dict[gid][other] = {'ani', 'af'}
            Hits of user genomes to GTDB representatives and to each other.
        """
        if self.qc_genomes is None:
            raise GTDBTkExit('parse() must be called before cluster().')

        # hit of ANI-screened genomes to their GTDB rep, so the final pass can compare it to de novo reps
        # (the GTDB hits of the other genomes are added by the caller)
        for gid, (gtdb_rep, hit) in gtdb_assignment.items():
            self.add_hit(ani_results, gid, gtdb_rep, hit)

        genomes_to_cluster = sorted(self.user_genomes - set(gtdb_assignment.keys()))
        self.logger.info(f'Dereplicating {len(genomes_to_cluster):,} genomes.')
        sp_cluster = dict()
        rounds = 0
        while True:
            rounds += 1
            sp_cluster = self.species_cluster(genomes_to_cluster, ani_results, sp_cluster)
            sp_cluster, gtdb_members, changing_reps_genomes, unassigned_genomes = \
                self.recluster_genomes(sp_cluster, gtdb_assignment, ani_results)
            # leftovers passing the quality gate go back to species_cluster;
            # each round adds at least one representative, so the loop ends
            genomes_to_cluster = sorted(gid for gid in unassigned_genomes
                                        if self.qc_genomes[gid].can_be_representative())
            if not genomes_to_cluster:
                break

        self.logger.info(f'{len(sp_cluster):,} species clusters found '
                         f'({rounds} round(s)).')
        if changing_reps_genomes:
            self.logger.info(f'{len(changing_reps_genomes):,} genome(s) assigned to a GTDB species '
                             f'by the ANI screen have been reclustered.')
        return sp_cluster, gtdb_members, changing_reps_genomes, unassigned_genomes

    @staticmethod
    def species_name(cluster):
        """ Name of a de novo species cluster."""
        return CONFIG.DEFAULT_SPECIES_CLUSTERS.format(inc=cluster['increment'])

    def write_clusters(self, path, sp_cluster, changing_reps_genomes):
        """ Write one line per user genome in a de novo species cluster."""
        with open(path, 'w') as fh:
            fh.write('user_genome\trepresentative\tde_novo_species\tradius\t'
                     'ani_to_representative\taf_to_representative\tprevious_gtdb_representative\n')
            for rep_id, cluster in sorted(sp_cluster.items(), key=lambda x: x[1]['increment']):
                species = self.species_name(cluster)
                previous = changing_reps_genomes.get(rep_id, ('N/A',))[0]
                fh.write(f"{rep_id}\t{rep_id}\t{species}\t{cluster['radius']}\tN/A\tN/A\t{previous}\n")
                for gid, hit in cluster['members']:
                    previous = changing_reps_genomes.get(gid, ('N/A',))[0]
                    fh.write(f"{gid}\t{rep_id}\t{species}\t{cluster['radius']}\t"
                             f"{round(hit['ani'], 2)}\t{round(hit['af'], 3)}\t{previous}\n")

    @staticmethod
    def gtdb_assignment_from_rows(list_summary_rows: List[ClassifySummaryFileRow]) -> Dict[str, Tuple[str, dict]]:
        """ GTDB assignment of the genomes classified by the ANI screen."""
        return {row.gid: (row.closest_genome_ref, {'ani': row.closest_genome_ani, 'af': row.closest_genome_af})
                for row in list_summary_rows}

    def __detect_format_checkm(self):
        with open(self.de_novo_file_path) as f:
            header = f.readline()
        sep = self.__detect_delimiters(header)
        both_models = {}
        header_infos = [h.strip() for h in header.strip().split(sep)]
        if all(item in header_infos for item in self.checkm2_headers_allmodels) and self.checkm2_headers_general in header_infos and self.checkm2_headers_specific in header_infos:
            self.logger.info('CheckM 2 file format detected')
            self.de_novo_file_format = 2
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
            raise GTDBTkExit('Unknown format of Completeness/Contamination file. the format should be either CheckM1, CheckM2 '
                             'or contains the headers "Name,Completeness and Contamination".')

    def __detect_delimiters(self,header):
        delimiters = ['\t', ',']
        for delimiter in delimiters:
            if delimiter in header:
                return delimiter
        return None
