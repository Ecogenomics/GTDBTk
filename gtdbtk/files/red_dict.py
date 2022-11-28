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
from gtdbtk.config.config import RED_DIST_ARC_DICT, RED_DIST_BAC_DICT
from gtdbtk.config.output import PATH_AR53_RED_DICT, PATH_BAC120_RED_DICT


class REDDictFile(object):
    """Store the GTDB-Tk RED dictionary."""

    def __init__(self, path: str, red_dict: Dict[str, float]):
        self.path = path
        self.data = red_dict

    def write(self):
        """Write the file to disk, note that domain is omitted."""
        make_sure_path_exists(os.path.dirname(self.path))
        with open(self.path, 'w') as fh:
            fh.write(f'Phylum\t{self.data["p__"]}\n')
            fh.write(f'Class\t{self.data["c__"]}\n')
            fh.write(f'Order\t{self.data["o__"]}\n')
            fh.write(f'Family\t{self.data["f__"]}\n')
            fh.write(f'Genus\t{self.data["g__"]}\n')


class REDDictFileAR53(REDDictFile):
    """Store the RED dictionary for the AR53 marker set."""

    def __init__(self, out_dir: str, prefix: str):
        path = os.path.join(out_dir, PATH_AR53_RED_DICT.format(prefix=prefix))
        super().__init__(path, RED_DIST_ARC_DICT)


class REDDictFileBAC120(REDDictFile):
    """Store the RED dictionary for the BAC120 marker set."""

    def __init__(self, out_dir: str, prefix: str):
        path = os.path.join(out_dir, PATH_BAC120_RED_DICT.format(prefix=prefix))
        super().__init__(path, RED_DIST_BAC_DICT)
