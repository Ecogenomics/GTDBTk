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

from gtdbtk.config.output import PATH_AR122_SUMMARY_OUT, PATH_BAC120_SUMMARY_OUT
from gtdbtk.exceptions import GTDBTkExit


class ClassifySummaryFileRow(object):
    """A row contained within the ClassifySummaryFile object."""

    def __init__(self):
        """Initialise the row, default all of the values to None."""
        self.gid = None
        self.classification = None
        self.fastani_ref = None
        self.fastani_ref_radius = None
        self.fastani_tax = None
        self.fastani_ani = None
        self.fastani_af = None
        self.closest_placement_ref = None
        self.closest_placement_radius = None
        self.closest_placement_tax = None
        self.closest_placement_ani = None
        self.closest_placement_af = None
        self.pplacer_tax = None
        self.classification_method = None
        self.note = None
        self.other_related_refs = None
        self.msa_percent = None
        self.tln_table = None
        self.red_value = None
        self.warnings = None


class ClassifySummaryFile(object):
    """Store the GTDB-Tk classify summary output."""

    def __init__(self, path: str, marker_set: str):
        self.logger = logging.getLogger('timestamp')
        self.path = path
        self.marker_set = marker_set
        self.rows = dict()  # keyed by user_genome
        self.none_value = 'N/A'

    @staticmethod
    def get_col_order(row: ClassifySummaryFileRow = None):
        """Return the column order that will be written. If a row is provided
        then format the row in that specific order."""
        if row is None:
            row = ClassifySummaryFileRow()
        mapping = [('user_genome', row.gid),
                   ('classification', row.classification),
                   ('fastani_reference', row.fastani_ref),
                   ('fastani_reference_radius', row.fastani_ref_radius),
                   ('fastani_taxonomy', row.fastani_tax),
                   ('fastani_ani', row.fastani_ani),
                   ('fastani_af', row.fastani_af),
                   ('closest_placement_reference', row.closest_placement_ref),
                   ('closest_placement_radius', row.closest_placement_radius),
                   ('closest_placement_taxonomy', row.closest_placement_tax),
                   ('closest_placement_ani', row.closest_placement_ani),
                   ('closest_placement_af', row.closest_placement_af),
                   ('pplacer_taxonomy', row.pplacer_tax),
                   ('classification_method', row.classification_method),
                   ('note', row.note),
                   ('other_related_references(genome_id,species_name,radius,ANI,AF)',
                    row.other_related_refs),
                   ('msa_percent', row.msa_percent),
                   ('translation_table', row.tln_table),
                   ('red_value', row.red_value),
                   ('warnings', row.warnings)]
        cols, data = list(), list()
        for col_name, col_val in mapping:
            cols.append(col_name)
            data.append(col_val)
        return cols, data

    def add_row(self, row: ClassifySummaryFileRow):
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
            raise GTDBTkExit(f'Error, classify summary file not found: {self.path}')
        with open(self.path) as fh:

            # Load and verify the columns match the expected order.
            cols_exp, _ = self.get_col_order()
            cols_cur = fh.readline().strip().split('\t')
            if cols_exp != cols_cur:
                raise GTDBTkExit(f'The classify summary file columns are inconsistent: {cols_cur}')

            # Process the data.
            for line in fh.readlines():
                data = line.strip().split('\t')
                row = ClassifySummaryFileRow()
                row.gid = data[0]
                row.classification = data[1]
                row.fastani_ref = data[2]
                row.fastani_ref_radius = data[3]
                row.fastani_tax = data[4]
                row.fastani_ani = data[5]
                row.fastani_af = data[6]
                row.closest_placement_ref = data[7]
                row.closest_placement_radius = data[8]
                row.closest_placement_tax = data[9]
                row.closest_placement_ani = data[10]
                row.closest_placement_af = data[11]
                row.pplacer_tax = data[12]
                row.classification_method = data[13]
                row.note = data[14]
                row.other_related_refs = data[15]
                row.msa_percent = data[16]
                row.tln_table = data[17]
                row.red_value = data[18]
                row.warnings = data[19]
                self.add_row(row)


class ClassifySummaryFileAR122(ClassifySummaryFile):
    """Store classify summary information for AR122 markers."""

    def __init__(self, out_dir: str, prefix: str):
        path = os.path.join(out_dir, PATH_AR122_SUMMARY_OUT.format(prefix=prefix))
        super().__init__(path, 'ar122')


class ClassifySummaryFileBAC120(ClassifySummaryFile):
    """Store classify summary information for BAC120 markers."""

    def __init__(self, out_dir: str, prefix: str):
        path = os.path.join(out_dir, PATH_BAC120_SUMMARY_OUT.format(prefix=prefix))
        super().__init__(path, 'bac120')
