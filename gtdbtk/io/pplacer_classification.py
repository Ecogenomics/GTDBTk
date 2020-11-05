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
from gtdbtk.config.output import PATH_AR122_PPLACER_CLASS, PATH_BAC120_PPLACER_CLASS
from gtdbtk.exceptions import GTDBTkExit


class PplacerClassifyFile(object):
    """Store the pplacer classifications."""

    def __init__(self, path: str):
        self.path = path
        self.data = dict()

    def add_genome(self, gid: str, tax_str: str):
        """Adds the pplacer classification of a given genome."""
        if gid in self.data:
            raise GTDBTkExit(f'Warning! Attempting to add duplicate genome: {gid}')
        self.data[gid] = tax_str

    def write(self):
        """Write the file to disk."""
        make_sure_path_exists(os.path.dirname(self.path))
        with open(self.path, 'w') as fh:
            for gid, tax_str in sorted(self.data.items()):
                fh.write(f'{gid}\t{tax_str}\n')


class PplacerClassifyFileAR122(PplacerClassifyFile):
    """Store the pplacer classifications for the AR122 marker set."""

    def __init__(self, out_dir: str, prefix: str):
        path = os.path.join(out_dir, PATH_AR122_PPLACER_CLASS.format(prefix=prefix))
        super().__init__(path)


class PplacerClassifyFileBAC120(PplacerClassifyFile):
    """Store the pplacer classifications for the BAC120 marker set."""

    def __init__(self, out_dir: str, prefix: str):
        path = os.path.join(out_dir, PATH_BAC120_PPLACER_CLASS.format(prefix=prefix))
        super().__init__(path)
