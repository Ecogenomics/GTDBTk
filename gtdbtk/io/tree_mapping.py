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

from gtdbtk.config.output import PATH_BAC120_TREE_MAPPING
from gtdbtk.exceptions import GTDBTkExit


class GenomeMappingFileRow(object):
    """A row contained within the GenomeMappingFile object."""

    def __init__(self):
        self.gid = None
        self.ani_classification = None
        self.mapped_tree = None
        self.rule = None

class GenomeMappingFile(object):
    """Store the GTDB-Tk classify summary output."""

    def __init__(self, out_dir: str,prefix: str):
        self.logger = logging.getLogger('timestamp')
        self.path = os.path.join(out_dir, PATH_BAC120_TREE_MAPPING.format(prefix=prefix))
        self.rows = dict()  # keyed by user_genome
        self.none_value = 'N/A'

    @staticmethod
    def get_col_order(row: GenomeMappingFileRow = None):
        """Return the column order that will be written. If a row is provided
        then format the row in that specific order."""
        if row is None:
            row = GenomeMappingFileRow()
        mapping = [('user_genome', row.gid),
                   ('is_ani_classification', row.ani_classification),
                   ('class_tree_mapped', row.mapped_tree),
                   ('classification_rule', row.rule)]
        cols, data = list(), list()
        for col_name, col_val in mapping:
            cols.append(col_name)
            data.append(col_val)
        return cols, data

    def add_row(self, row: GenomeMappingFileRow):
        if row.gid in self.rows:
            raise GTDBTkExit(f'Attempting to add duplicate row: {row.gid}')
        self.rows[row.gid] = row

    def write(self):
        """Writes the summary file to disk. None will be replaced with N/A"""
        with open(self.path, 'w') as fh:
            fh.write('\t'.join(self.get_col_order()[0]) + '\n')
            for gid, row in sorted(self.rows.items()):
                buf = list()
                for data in self.get_col_order(row)[1]:
                    buf.append(self.none_value if data is None else str(data))
                fh.write('\t'.join(buf) + '\n')

    def read(self):
        """Read the summary file from disk."""
        if not os.path.isfile(self.path):
            raise GTDBTkExit(f'Error, classify tree mappings file not found: {self.path}')
        with open(self.path) as fh:

            # Load and verify the columns match the expected order.
            cols_exp, _ = self.get_col_order()
            cols_cur = fh.readline().strip().split('\t')
            if cols_exp != cols_cur:
                raise GTDBTkExit(f'The classify tree mappings columns are inconsistent: {cols_cur}')

            # Process the data.
            for line in fh.readlines():
                data = line.strip().split('\t')
                row = GenomeMappingFileRow()
                row.gid = data[0]
                row.ani_classification = data[1]
                row.mapped_tree = data[2]
                row.rule = data[3]
                self.add_row(row)