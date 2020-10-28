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

from gtdbtk.biolib_lite.common import make_sure_path_exists
from gtdbtk.config.output import PATH_TLN_TABLE_SUMMARY
from gtdbtk.exceptions import GTDBTkExit


class TlnTableSummaryFile(object):
    """Records the translation table for one or more genomes."""

    def __init__(self, out_dir: str, prefix: str):
        """Configure paths and initialise storage dictionary."""
        self.path = os.path.join(out_dir, PATH_TLN_TABLE_SUMMARY.format(prefix=prefix))
        self.genomes = dict()

    def add_genome(self, genome_id: str, tln_table: int):
        """Record a translation table for a genome."""
        if genome_id in self.genomes:
            raise GTDBTkExit(f'Genome already exists in summary file: {genome_id}')
        self.genomes[genome_id] = tln_table

    def write(self):
        """Write the translation table summary file to disk."""
        make_sure_path_exists(os.path.dirname(self.path))
        with open(self.path, 'w') as fh:
            for genome_id, tln_table in sorted(self.genomes.items()):
                fh.write(f'{genome_id}\t{tln_table}\n')

    def read(self):
        """Read the translation table summary file from disk."""
        if len(self.genomes) > 0:
            raise GTDBTkExit(f'Warning! Attempting to override in-memory values '
                             f'for translation table summary file: {self.path}')
        with open(self.path, 'r') as fh:
            for line in fh.readlines():
                gid, tbl = line.strip().split('\t')
                self.genomes[gid] = int(tbl)
