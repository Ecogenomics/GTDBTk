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
from typing import Dict

from gtdbtk.biolib_lite.common import make_sure_path_exists
from gtdbtk.config.output import PATH_AR53_DISAPPEARING_GENOMES, PATH_BAC120_DISAPPEARING_GENOMES
from gtdbtk.exceptions import GTDBTkExit


class DisappearingGenomesFile(object):
    """Store the GTDB-Tk RED dictionary."""

    def __init__(self, path: str,domain: str):
        self.path: str = path
        self.domain: str = domain
        self.data: Dict[str, str] = dict()
        self.file_name: str = os.path.basename(path)

    def add_genome(self, gid: str, tree_index: str):
        """PlAdds the pplacer classification of a given genome."""
        if gid in self.data:
            raise GTDBTkExit(f'Warning! Attempting to add duplicate genome: {gid}')
        self.data[gid] = tree_index

    def write(self):
        """Write the file to disk, note that domain is omitted."""
        make_sure_path_exists(os.path.dirname(self.path))
        with open(self.path, 'w') as fh:
            fh.write(f'genome_id\tdomain\ttree_index\n')
            for seqid, infos in self.data.items():
                fh.write(f'{seqid}\t{self.domain}\t{infos}\n')


class DisappearingGenomesFileAR53(DisappearingGenomesFile):
    """Store the RED dictionary for the AR53 marker set."""

    def __init__(self, out_dir: str, prefix: str):
        path = os.path.join(out_dir, PATH_AR53_DISAPPEARING_GENOMES.format(prefix=prefix))
        super().__init__(path,'archaea')


class DisappearingGenomesFileBAC120(DisappearingGenomesFile):
    """Store the RED dictionary for the BAC120 marker set."""

    def __init__(self, out_dir: str, prefix: str):
        path = os.path.join(out_dir, PATH_BAC120_DISAPPEARING_GENOMES.format(prefix=prefix))
        super().__init__(path,'bacteria')
