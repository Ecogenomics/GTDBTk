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
from gtdbtk.config.output import PATH_AR53_PPLACER_CLASS, PATH_BAC120_PPLACER_CLASS, PATH_BAC120_BACKBONE_PPLACER_CLASS, \
    PATH_BAC120_CLASS_LEVEL_PPLACER_CLASS
from gtdbtk.exceptions import GTDBTkExit


class PplacerClassifyFile(object):
    """Store the pplacer classifications."""

    def __init__(self, path: str):
        self.path: str = path
        self.data: Dict[str, str] = dict()


    def add_genome(self, gid: str, tax_str: str):
        """Adds the pplacer classification of a given genome."""
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


class PplacerClassifyFileAR53(PplacerClassifyFile):
    """Store the pplacer classifications for the AR53 marker set."""

    def __init__(self, out_dir: str, prefix: str):
        path = os.path.join(out_dir, PATH_AR53_PPLACER_CLASS.format(prefix=prefix))
        super().__init__(path)


class PplacerClassifyFileBAC120(PplacerClassifyFile):
    """Store the pplacer classifications for the BAC120 marker set."""

    def __init__(self, out_dir: str, prefix: str):
        path = os.path.join(out_dir, PATH_BAC120_PPLACER_CLASS.format(prefix=prefix))
        super().__init__(path)

class PplacerLowClassifyFileBAC120(PplacerClassifyFile):
    """Store the pplacer classifications for the BAC120 marker set."""

    def __init__(self, out_dir: str, prefix: str,iter:str):
        path = os.path.join(out_dir, PATH_BAC120_CLASS_LEVEL_PPLACER_CLASS.format(prefix=prefix,iter=iter))
        super().__init__(path)


class PplacerHighClassifyRow(object):
    """Initialise the row, default all of the values to None."""

    def __init__(self):
        self.gid = None
        self.gtdb_taxonomy_red = None
        self.gtdb_taxonomy_terminal = None
        self.pplacer_taxonomy = None
        self.is_terminal = None
        self.red = None


class PplacerHighClassifyFile(object):
    """Store the pplacer classifications."""

    def __init__(self,out_dir: str,prefix: str):
        self.path = os.path.join(out_dir, PATH_BAC120_BACKBONE_PPLACER_CLASS.format(prefix=prefix))
        self.rows = dict()  # keyed by user_genome
        self.none_value = 'N/A'

    def add_row(self, row: PplacerHighClassifyRow):
        if row.gid in self.rows:
            raise GTDBTkExit(f'Attempting to add duplicate row: {row.gid}')
        self.rows[row.gid] = row

    @staticmethod
    def get_col_order(row: PplacerHighClassifyRow = None):
        """Return the column order that will be written. If a row is provided
        then format the row in that specific order."""
        if row is None:
            row = PplacerHighClassifyRow()
        mapping = [('user_genome', row.gid),
                   ('gtdb_taxonomy_red', row.gtdb_taxonomy_red),
                   ('gtdb_taxonomy_terminal', row.gtdb_taxonomy_terminal),
                   ('pplacer_taxonomy', row.pplacer_taxonomy),
                   ('is_terminal', row.is_terminal),
                   ('red', row.red)]
        cols, data = list(), list()
        for col_name, col_val in mapping:
            cols.append(col_name)
            data.append(col_val)
        return cols, data


    def write(self):
        """Write the file to disk."""
        make_sure_path_exists(os.path.dirname(self.path))
        cols = ['gid','gtdb_taxonomy','pplacer_taxonomy','is_terminal','red']
        with open(self.path, 'w') as fh:
            fh.write('\t'.join(self.get_col_order()[0]) + '\n')
            for gid, row in sorted(self.rows.items()):
                buf = list()
                for data in self.get_col_order(row)[1]:
                    buf.append(self.none_value if data is None else str(data))
                fh.write('\t'.join(buf) + '\n')
