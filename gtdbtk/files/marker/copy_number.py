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
from typing import Set, Dict, List, Union

from gtdbtk.biolib_lite.common import make_sure_path_exists
from gtdbtk.biolib_lite.seq_io import read_fasta
from gtdbtk.config.config import BAC120_MARKERS, AR53_MARKERS
from gtdbtk.config.output import PATH_BAC120_MARKER_SUMMARY, PATH_AR53_MARKER_SUMMARY
from gtdbtk.exceptions import GTDBTkExit
from gtdbtk.files.marker.tophit import TopHitPfamFile, TopHitTigrFile, Hit


class CopyNumberFile(object):
    """Store copy number information for a specific marker set. Records the
    sequence of the (most significant) hits for the alignment step."""

    def __init__(self, path: str, marker_set: str, marker_dict: Dict[str, List[str]]):
        """Initialise the file and extract all of the marker names."""
        self.logger = logging.getLogger('timestamp')
        self.path = path
        self.marker_set = marker_set
        self.genomes = dict()
        self.marker_names = self._extract_marker_names(marker_dict)

    def add_genome(self, genome_id: str, path_faa: str, pfam_th: TopHitPfamFile, tigr_th: TopHitTigrFile):
        """Process the top hit files for a genome and store the copy info."""
        if genome_id in self.genomes:
            self.logger.warning(f'Genome already exists in copy number file: {genome_id}')
        self.genomes[genome_id] = {'unq': dict(), 'mul': dict(), 'muq': dict(), 'mis': dict()}

        # Pointers to unique, multiple hit, multiple-unique, missing markers.
        cur_unq = self.genomes[genome_id]['unq']
        cur_mul = self.genomes[genome_id]['mul']
        cur_muq = self.genomes[genome_id]['muq']
        cur_mis = self.genomes[genome_id]['mis']

        # Load genes from the prodigal faa file.
        d_genes = read_fasta(path_faa, False)
        for seq_id, seq in d_genes.items():
            if seq.endswith('*'):
                d_genes[seq_id] = seq[:-1]

        # Create a dictionary of marker names -> Hits
        d_hmm_hits = self._merge_hit_files(pfam_th, tigr_th)

        # Foreach expected marker determine which category it falls into.
        for marker_id in self.marker_names:

            # Marker is missing.
            if marker_id not in d_hmm_hits:
                cur_mis[marker_id] = None

            # Multiple hits to to the same marker.
            elif len(d_hmm_hits[marker_id]) > 1:

                # If sequences are the same, take the most significant hit
                unq_seqs = {d_genes[x.gene_id] for x in d_hmm_hits[marker_id]}
                if len(unq_seqs) == 1:
                    cur_top_hit = sorted(d_hmm_hits[marker_id], reverse=True)[0]
                    cur_muq[marker_id] = {'hit': cur_top_hit, 'seq': d_genes[cur_top_hit.gene_id]}

                # Marker maps to multiple genes.
                else:
                    cur_mul[marker_id] = None

            # This was a unique hit.
            else:
                cur_hit = d_hmm_hits[marker_id][0]
                cur_unq[marker_id] = {'hit': cur_hit, 'seq': d_genes[cur_hit.gene_id]}

        # Sanity check - confirm that the total number of markers matches.
        if len(self.marker_names) != len(cur_unq) + len(cur_mul) + len(cur_muq) + len(cur_mis):
            raise GTDBTkExit('The marker set is inconsistent, please report this issue.')

    @staticmethod
    def _extract_marker_names(marker_dict: Dict[str, List[str]]) -> Set[str]:
        """Parse the GTDB-Tk configuration file to get the HMM names."""

        def __remove_suffix(in_str: str) -> str:
            return in_str.replace('.HMM', '').replace('.hmm', '')

        out = set()
        for m_list in marker_dict.values():
            for m_name in m_list:
                out.add(__remove_suffix(m_name))
        return out

    @staticmethod
    def _merge_hit_files(pfam_th: TopHitPfamFile, tigr_th: TopHitTigrFile) -> Dict[str, List[Hit]]:
        """Returns a list of Hits grouped by HMM id as a dictionary."""
        out = dict()
        for cur_tophit in (pfam_th, tigr_th):
            for gene_id, cur_hit in cur_tophit.iter_hits():
                if cur_hit.hmm_id not in out:
                    out[cur_hit.hmm_id] = list()
                out[cur_hit.hmm_id].append(cur_hit)
        return out

    def get_single_copy_hits(self, genome_id: str) -> Dict[str, Dict[str, Union[str, Hit]]]:
        """Return hit and seq info for single copy markers."""
        return {**self.genomes[genome_id]['unq'], **self.genomes[genome_id]['muq']}

    def write(self):
        """Write the summary file to disk."""
        make_sure_path_exists(os.path.dirname(self.path))
        header = ['name', 'number_unique_genes', 'number_multiple_genes',
                  'number_multiple_unique_genes', 'number_missing_genes',
                  'list_unique_genes', 'list_multiple_genes',
                  'list_multiple_unique_genes', 'list_missing_genes']
        with open(self.path, 'w') as fh:
            fh.write('\t'.join(header) + '\n')
            for genome_id, marker_dict in sorted(self.genomes.items()):
                fh.write(f'{genome_id}\t'
                         f'{len(marker_dict["unq"])}\t'
                         f'{len(marker_dict["mul"])}\t'
                         f'{len(marker_dict["muq"])}\t'
                         f'{len(marker_dict["mis"])}\t'
                         f'{",".join(sorted(marker_dict["unq"]))}\t'
                         f'{",".join(sorted(marker_dict["mul"]))}\t'
                         f'{",".join(sorted(marker_dict["muq"]))}\t'
                         f'{",".join(sorted(marker_dict["mis"]))}\n')

    def read(self):
        """Reads the marker names from disk. No sequence information!"""
        with open(self.path) as fh:
            fh.readline()
            for line in fh.readlines():
                genome_id, n_unq, n_mul, n_muq, n_mis, unq, mul, muq, mis = line.split('\t')
                n_unq, n_mul, n_muq, n_mis = int(n_unq), int(n_mul), int(n_muq), int(n_mis)
                self.genomes[genome_id] = {'unq': {x: None for x in unq.strip().split(',') if len(x) > 0},
                                           'mul': {x: None for x in mul.strip().split(',') if len(x) > 0},
                                           'muq': {x: None for x in muq.strip().split(',') if len(x) > 0},
                                           'mis': {x: None for x in mis.strip().split(',') if len(x) > 0}}
                cur_dict = self.genomes[genome_id]
                if len(cur_dict['unq']) != n_unq or len(cur_dict['mul']) != n_mul or \
                        len(cur_dict['muq']) != n_muq or len(cur_dict['mis']) != n_mis:
                    raise GTDBTkExit(f'The marker file is inconsistent: {self.path}')


class CopyNumberFileAR53(CopyNumberFile):
    """Store hmm copy number information for AR53 markers."""

    def __init__(self, out_dir: str, prefix: str):
        path = os.path.join(out_dir, PATH_AR53_MARKER_SUMMARY.format(prefix=prefix))
        super().__init__(path, 'ar53', AR53_MARKERS)


class CopyNumberFileBAC120(CopyNumberFile):
    """Store hmm copy number information for BAC120 markers."""

    def __init__(self, out_dir: str, prefix: str):
        path = os.path.join(out_dir, PATH_BAC120_MARKER_SUMMARY.format(prefix=prefix))
        super().__init__(path, 'bac120', BAC120_MARKERS)
