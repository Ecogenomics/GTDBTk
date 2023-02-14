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
from typing import Dict, List, Tuple, Optional, Union

from gtdbtk.biolib_lite.common import make_sure_path_exists
from gtdbtk.config.output import PATH_AR53_SUMMARY_OUT, PATH_BAC120_SUMMARY_OUT
from gtdbtk.exceptions import GTDBTkExit


class ClassifySummaryFileRow:
    """A row contained within the ClassifySummaryFile object."""

    __slots__ = ('gid', 'classification', 'fastani_ref', 'fastani_ref_radius',
                 'fastani_tax', 'fastani_ani', 'fastani_af', 'closest_placement_ref',
                 'closest_placement_radius', 'closest_placement_tax', 'closest_placement_ani',
                 'closest_placement_af', 'pplacer_tax', 'classification_method',
                 'note', 'other_related_refs', 'msa_percent', 'tln_table',
                 'red_value', 'warnings')

    def __init__(self):
        """Initialise the row, default all the values to None."""
        self.gid: Optional[str] = None
        self.classification: Optional[str] = None
        self.fastani_ref: Optional[str] = None
        self.fastani_ref_radius: Optional[float] = None
        self.fastani_tax: Optional[str] = None
        self.fastani_ani: Optional[float] = None
        self.fastani_af: Optional[float] = None
        self.closest_placement_ref: Optional[str] = None
        self.closest_placement_radius: Optional[float] = None
        self.closest_placement_tax: Optional[str] = None
        self.closest_placement_ani: Optional[float] = None
        self.closest_placement_af: Optional[float] = None
        self.pplacer_tax: Optional[str] = None
        self.classification_method: Optional[str] = None
        self.note: Optional[str] = None
        self.other_related_refs: Optional[str] = None
        self.msa_percent: Optional[float] = None
        self.tln_table: Optional[int] = None
        self.red_value: Optional[float] = None
        self.warnings: Optional[str] = None


class ClassifySummaryFile:
    """Store the GTDB-Tk classify summary output."""

    __slots__ = ('logger', 'path', 'marker_set', 'rows', 'none_value')

    def __init__(self, path: str, marker_set: Optional[str] = None):
        self.logger = logging.getLogger('timestamp')
        self.path: str = path
        self.marker_set: Optional[str] = marker_set
        self.rows: Dict[str, ClassifySummaryFileRow] = dict()  # keyed by user_genome
        self.none_value: str = 'N/A'

    @staticmethod
    def get_col_order(row: ClassifySummaryFileRow = None) -> Tuple[List[str], List[Union[str, float, int]]]:
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
                   ('other_related_references(genome_id,species_name,radius,ANI,AF)', row.other_related_refs),
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

    def has_row(self) -> bool:
        if self.rows.items():
            return True
        return False

    def get_gid_taxonomy(self) -> Dict[str, List[str]]:
        out = dict()
        for gid, row in self.rows.items():
            split_tax = row.classification.split(';')
            if len(split_tax) != 7:
                raise GTDBTkExit(f'Expected a 7-rank taxonomy for {gid} but got {row.classification}')
            out[gid] = split_tax
        return out

    def write(self):
        """Writes the summary file to disk. None will be replaced with N/A"""
        make_sure_path_exists(os.path.dirname(self.path))
        with open(self.path, 'w') as fh:
            fh.write('\t'.join(self.get_col_order()[0]) + '\n')
            for gid, row in sorted(self.rows.items()):
                buf = list()
                for idx,data in enumerate(self.get_col_order(row)[1]):
                    # for the red_value field, we want to round the data to 5 decimals after the comma if the value is not None
                    if idx==self.get_col_order()[0].index('red_value') and data is not None:
                        data = round(data,5)
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


class ClassifySummaryFileAR53(ClassifySummaryFile):
    """Store classify summary information for AR53 markers."""

    def __init__(self, out_dir: str, prefix: str):
        path = os.path.join(out_dir, PATH_AR53_SUMMARY_OUT.format(prefix=prefix))
        super().__init__(path, 'ar53')


class ClassifySummaryFileBAC120(ClassifySummaryFile):
    """Store classify summary information for BAC120 markers."""

    def __init__(self, out_dir: str, prefix: str):
        path = os.path.join(out_dir, PATH_BAC120_SUMMARY_OUT.format(prefix=prefix))
        super().__init__(path, 'bac120')
