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
from shutil import copyfile

import gtdbtk.config.config as Config
from gtdbtk.biolib_lite.common import make_sure_path_exists
from gtdbtk.biolib_lite.execute import check_dependencies
from gtdbtk.biolib_lite.logger import colour
from gtdbtk.biolib_lite.seq_io import read_fasta
from gtdbtk.exceptions import GTDBTkException, GTDBTkExit
from gtdbtk.tools import sha1_dir


class Misc(object):

    def __init__(self):
        """Initialize."""
        self.logger = logging.getLogger('timestamp')

    def trim_msa(self, untrimmed_msa, mask_type, maskid, output_file):
        """Trim the multiple sequence alignment using a mask.

        Parameters
        ----------
        untrimmed_msa : str
            The path to the untrimmed MSA.
        mask_type : str
            Which mask should be used, reference or user specified.
        maskid : str
            The path to the mask used for trimming.
        output_file : str
            The path to the output trimmed MSA.
        """
        if maskid == 'bac' and mask_type == 'reference':
            mask = os.path.join(Config.MASK_DIR, Config.MASK_BAC120)
        elif maskid == 'arc' and mask_type == 'reference':
            mask = os.path.join(Config.MASK_DIR, Config.MASK_AR122)
        elif mask_type == 'file':
            mask = maskid
        else:
            self.logger.error('Command not understood.')
            raise GTDBTkException('Command not understood.')

        with open(mask, 'r') as f:
            maskstr = f.readline()

        with open(output_file, 'w') as outfwriter:
            dict_genomes = read_fasta(untrimmed_msa, False)

            for k, v in dict_genomes.items():
                aligned_seq = ''.join([v[i] for i in range(0, len(maskstr)) if maskstr[i] == '1'])
                fasta_outstr = ">%s\n%s\n" % (k, aligned_seq)
                outfwriter.write(fasta_outstr)

    def export_msa(self, domain, output_file):
        """Export the MSA to a file, create the path if it doesn't exist.

        Parameters
        ----------
        domain : str
            The domain used to determine the marker set.
        output_file : str
            The path where the MSA should be exported.
        """
        file_to_export = Config.CONCAT_BAC120
        if domain == 'arc':
            file_to_export = Config.CONCAT_AR122

        make_sure_path_exists(os.path.dirname(output_file))
        copyfile(file_to_export, output_file)

    def checkfile(self, file_path, file_name):
        """Check that a file exists, output the result to the logger.

        Returns
        -------
        bool
            True if the folder exists, False otherwise.
        """
        if os.path.exists(file_path) and os.path.getsize(file_path) > 0:
            self.logger.info("Check file {} ({}): {}".format(
                file_name, file_path, colour('OK', ['bright'], fg='green')))
            return True
        else:
            self.logger.warning("Check file {} ({}): {}".format(
                file_name, file_path, colour('MISSING', ['bright'], fg='red')))
            return False

    def checkfolder(self, folder_path, folder_name):
        """Check that a folder exists, output the result to the logger.

        Returns
        -------
        bool
            True if the folder exists, False otherwise.
        """
        if os.path.isdir(folder_path) and len(os.listdir(folder_path)) > 0:
            self.logger.info("Check folder {} ({}): {}".format(
                folder_name, folder_path, colour('OK', ['bright'], fg='green')))
            return True
        else:
            self.logger.warning("Check folder {} ({}): {}".format(
                folder_name, folder_path, colour('MISSING', ['bright'], fg='red')))
            return False

    def check_install(self):
        """Check that all reference files exist.

        Returns
        -------
        bool
            True if the installation is complete, False otherwise.
        """

        # Check that all programs are on the system path.
        self.logger.info(f'Checking that all third-party software are on the system path:')
        names = {'prodigal', 'hmmsearch', 'fastANI', 'mash', 'pplacer', 'guppy',
                 'FastTree', 'FastTreeMP', 'hmmalign'}
        for name in sorted(names):
            on_path = False
            try:
                on_path = on_path or check_dependencies([name], exit_on_fail=False)
            except:
                pass
            if on_path:
                self.logger.info("         |-- {:16} {}".format(
                    name, colour('OK', ['bright'], fg='green')))
            else:
                self.logger.info("         |-- {:16} {}".format(
                    name, colour('NOT FOUND', ['bright'], fg='yellow')))

        # Assume this was successful unless otherwise observed.
        ok = True

        # Compute the hash for each directory
        self.logger.info(f'Checking integrity of reference package: {Config.GENERIC_PATH}')
        for obj_path, expected_hash in Config.REF_HASHES.items():
            base_name = obj_path[:-1] if obj_path.endswith('/') else obj_path
            base_name = base_name.split('/')[-1]
            user_hash = sha1_dir(obj_path, progress=True)

            if user_hash != expected_hash:
                self.logger.info("         |-- {:16} {}".format(
                    base_name, colour('HASH MISMATCH', ['bright'], fg='yellow')))
                ok = False
            else:
                self.logger.info("         |-- {:16} {}".format(
                    base_name, colour('OK', ['bright'], fg='green')))

        if not ok:
            raise GTDBTkExit('Unexpected files were seen, or the reference package is corrupt.')
