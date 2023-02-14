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
from typing import Tuple, Optional, Iterator

from gtdbtk.biolib_lite.common import make_sure_path_exists
from gtdbtk.config.output import TIGRFAM_TOP_HIT_SUFFIX, PFAM_TOP_HIT_SUFFIX, CHECKSUM_SUFFIX
from gtdbtk.tools import sha256


class Hit(object):
    """More significant hits have a higher bit-score and lower e-value.

    Running sorted() on a list of Hit objects will yield a list in order of
    least significant to most significant hits.
    """

    def __init__(self, gene_id: str, hmm_id: str, e_val: float, bit_score: float):
        """Store information about this hit."""
        self.gene_id = gene_id
        self.hmm_id = hmm_id
        self.e_val = e_val
        self.bit_score = bit_score

    def __repr__(self) -> str:
        return f'{self.gene_id} {self.hmm_id} ({self.e_val}/{self.bit_score})'

    def __eq__(self, other) -> bool:
        """Returns True if self is equal to other."""
        return isinstance(other, Hit) and other.gene_id == self.gene_id \
               and other.hmm_id == self.hmm_id and other.e_val == self.e_val \
               and other.bit_score == self.bit_score

    def __lt__(self, other) -> bool:
        """Is self less significant than other?"""
        if self.bit_score < other.bit_score:
            return True
        elif self.bit_score == other.bit_score:
            if self.e_val > other.e_val:
                return True
            elif self.e_val == other.e_val:
                if self.hmm_id > other.hmm_id:
                    return True
                elif self.hmm_id == other.hmm_id:
                    return self.gene_id > other.gene_id
        return False

    def __gt__(self, other) -> bool:
        """Is self more significant than other?"""
        raise NotImplemented
        # if self.bit_score > other.bit_score:
        #     return True
        # elif self.bit_score == other.bit_score:
        #     if self.e_val < other.e_val:
        #         return True
        #     elif self.e_val == other.e_val:
        #         if self.hmm_id < other.hmm_id:
        #             return True
        #         elif self.hmm_id == other.hmm_id:
        #             return self.gene_id < other.gene_id
        # return False

    def __hash__(self) -> int:
        return hash(f'{self.gene_id}_{self.hmm_id}_{self.e_val}_{self.bit_score}')

    def hmm_str(self) -> str:
        """Report the e-value and bit-score for the hmm."""
        return f'{self.hmm_id},{self.e_val},{self.bit_score}'


class TopHitFile(object):
    """Stores information about the top marker hits for Pfam/TIGRFAM markers."""

    def __init__(self, path: str):
        """Setup the path to the file and initialise storage dictionary."""
        self.path = path
        self.hits = dict()

    def add_hit(self, gene_id: str, hmm_id: str, e_val: float, bit_score: float):
        """Store the most significant HMM hit for each gene."""
        if gene_id not in self.hits:
            self.hits[gene_id] = dict()

        # Check if this hit already exists.
        new_hit = Hit(gene_id, hmm_id, e_val, bit_score)
        if hmm_id in self.hits[gene_id]:
            # Store the new hit if it's more significant.
            if self.hits[gene_id][hmm_id] < new_hit:
                self.hits[gene_id][hmm_id] = new_hit
        else:
            self.hits[gene_id][hmm_id] = new_hit

    def contains_gene_id(self, gene_id: str) -> bool:
        """Returns True if the gene_id is found in the top hit file."""
        return gene_id in self.hits

    def contains_gene_hmm(self, gene_id: str, hmm_id: str) -> bool:
        """Returns True if there is a hmm hit for the gene_id."""
        return gene_id in self.hits and hmm_id in self.hits[gene_id]

    def get_top_hit(self, gene_id: str) -> Optional[Hit]:
        """Returns the most significant hit for a given gene id or None."""
        if gene_id not in self.hits or len(self.hits[gene_id]) == 0:
            return None
        return sorted(self.hits[gene_id].values(), reverse=True)[0]

    def get_hmm_hit(self, gene_id: str, hmm_id: str) -> Hit:
        """Returns the hit indexed by the gene and hmm id."""
        return self.hits[gene_id][hmm_id]

    def write(self):
        """Writes the file to disk and creates a checksum."""
        # Write the top hit file.
        make_sure_path_exists(os.path.dirname(self.path))
        header = ['Gene Id', 'Top hits (Family id,e-value,bitscore)']
        with open(self.path, 'w') as fh:
            fh.write('\t'.join(header) + '\n')
            for gene_id, hits in sorted(self.hits.items()):
                out_hits = list()
                for cur_hit in sorted(hits.values(), reverse=True):
                    out_hits.append(cur_hit.hmm_str())
                concat_hits = ';'.join(out_hits)
                fh.write(f'{gene_id}\t{concat_hits}\n')

        # Write the checksum.
        with open(f'{self.path}{CHECKSUM_SUFFIX}', 'w') as fh:
            fh.write(sha256(self.path))

    def read(self):
        """Read the contents of an existing tophit file."""
        with open(self.path) as fh:
            fh.readline()
            for line in fh:
                gene_id, hits = line.strip().split('\t')
                for hit in hits.split(';'):
                    hmm_id, e_val, bit_score = hit.split(',')
                    e_val, bit_score = float(e_val), float(bit_score)
                    self.add_hit(gene_id, hmm_id, e_val, bit_score)

    def iter_hits(self) -> Iterator[Tuple[str, Hit]]:
        """Iterate over all genes and their respective hits."""
        for gene_id, hmm_dict in self.hits.items():
            for cur_hit in hmm_dict.values():
                yield gene_id, cur_hit


class TopHitPfamFile(TopHitFile):
    """The top hit file given the Pfam output."""

    def __init__(self, out_dir: str, gid: str):
        path = self.get_path(out_dir, gid)
        super().__init__(path)

    @staticmethod
    def get_path(out_dir: str, gid: str):
        return os.path.join(out_dir, gid, f'{gid}{PFAM_TOP_HIT_SUFFIX}')


class TopHitTigrFile(TopHitFile):
    """The top hit file given the Tigrfam output."""

    def __init__(self, out_dir: str, gid: str):
        path = self.get_path(out_dir, gid)
        super().__init__(path)

    @staticmethod
    def get_path(out_dir: str, gid: str):
        return os.path.join(out_dir, gid, f'{gid}{TIGRFAM_TOP_HIT_SUFFIX}')

    def add_hit(self, gene_id: str, hmm_id: str, e_val: float, bit_score: float):
        """Only store the most significant HMM hit per gene as per GTDBNCBI."""
        if self.contains_gene_id(gene_id):
            # The new hit was more significant, remove any other hits.
            if self.get_top_hit(gene_id) < Hit(gene_id, hmm_id, e_val, bit_score):
                self.hits[gene_id] = dict()
                super().add_hit(gene_id, hmm_id, e_val, bit_score)
        else:
            super().add_hit(gene_id, hmm_id, e_val, bit_score)
