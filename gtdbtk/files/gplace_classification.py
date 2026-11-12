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
from gtdbtk.config.output import PATH_AR53_GPLACE_CLASS, PATH_BAC120_GPLACE_CLASS
from gtdbtk.exceptions import GTDBTkExit


class GplaceClassifyFile(object):
    """Store the gPlace classifications."""

    def __init__(self, path: str):
        self.path: str = path
        self.data: Dict[str, str] = dict()


    def add_genome(self, gid: str, tax_str: str):
        """Adds the gplace classification of a given genome."""
        if gid in self.data:
            raise GTDBTkExit(f'Warning! Attempting to add duplicate genome: {gid}')
        self.data[gid] = tax_str

    def write(self):
        """Write the file to disk."""
        make_sure_path_exists(os.path.dirname(self.path))
        if len(self.data) > 0 :
            with open(self.path, 'w') as fh:
                for gid, tax_str in sorted(self.data.items()):
                    fh.write(f'{gid}\t{tax_str}\n')


class GplaceClassifyFileAR53(GplaceClassifyFile):
    """Store the gPlace classifications for the AR53 marker set."""

    def __init__(self, out_dir: str, prefix: str):
        path = os.path.join(out_dir, PATH_AR53_GPLACE_CLASS.format(prefix=prefix))
        super().__init__(path)



class GplaceClassifyFileBAC120(GplaceClassifyFile):
    """Store the gPlace classifications for the BAC120 marker set."""

    def __init__(self, out_dir: str, prefix: str):
        path = os.path.join(out_dir, PATH_BAC120_GPLACE_CLASS.format(prefix=prefix))
        super().__init__(path)





