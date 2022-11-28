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
from gtdbtk.config.config import AR53_MARKERS, BAC120_MARKERS, TIGRFAM_HMMS, PFAM_HMM_DIR
from gtdbtk.config.output import PATH_AR53_MARKER_INFO, PATH_BAC120_MARKER_INFO


class MarkerInfoFile(object):
    """Store the GTDB-Tk RED dictionary."""

    marker_paths = {"PFAM": os.path.join(PFAM_HMM_DIR, 'individual_hmms'),
                    "TIGRFAM": os.path.join(os.path.dirname(TIGRFAM_HMMS), 'individual_hmms')}

    def __init__(self, path: str, markers: dict):
        self.path = path
        self.markers = self._parse_markers(markers)

    def _parse_markers(self, markers):
        out = dict()
        for db_marker in sorted(markers):
            for marker in markers[db_marker]:
                marker_id = marker[0:marker.rfind('.')]
                marker_path = os.path.join(self.marker_paths[db_marker], marker)

                # get marker name, description, and size
                with open(marker_path) as fh:
                    for line in fh:
                        if line.startswith("NAME  "):
                            marker_name = line.split("  ")[1].strip()
                        elif line.startswith("DESC  "):
                            marker_desc = line.split("  ")[1].strip()
                        elif line.startswith("LENG  "):
                            marker_size = int(line.split("  ")[1].strip())
                            break

                out[marker_id] = {'name': marker_name,
                                  'desc': marker_desc,
                                  'size': marker_size,
                                  'path': marker_path}
        return out

    def write(self):
        """Write the file to disk, note that domain is omitted."""
        make_sure_path_exists(os.path.dirname(self.path))
        with open(self.path, 'w') as fh:
            fh.write('Marker ID\tName\tDescription\tLength (bp)\n')
            for marker_id, marker_d in sorted(self.markers.items()):
                row = [marker_id, marker_d['name'], marker_d['desc'], str(marker_d['size'])]
                fh.write('\t'.join(row) + '\n')


class MarkerInfoFileAR53(MarkerInfoFile):
    """Marker information for the AR53 marker set."""

    def __init__(self, out_dir: str, prefix: str):
        path = os.path.join(out_dir, PATH_AR53_MARKER_INFO.format(prefix=prefix))
        super().__init__(path, AR53_MARKERS)


class MarkerInfoFileBAC120(MarkerInfoFile):
    """Marker information for the BAC120 marker set."""

    def __init__(self, out_dir: str, prefix: str):
        path = os.path.join(out_dir, PATH_BAC120_MARKER_INFO.format(prefix=prefix))
        super().__init__(path, BAC120_MARKERS)
