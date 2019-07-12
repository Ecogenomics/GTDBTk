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

import config.config as Config
from biolib_lite.common import make_sure_path_exists
from biolib_lite.seq_io import read_fasta
from gtdbtk.exceptions import ReferenceFileMalformed


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
        with open(mask, 'r') as f:
            maskstr = f.readline()

        outfwriter = open(output_file, 'w')
        dict_genomes = read_fasta(untrimmed_msa, False)

        for k, v in dict_genomes.iteritems():
            aligned_seq = ''.join([v[i] for i in range(0, len(maskstr)) if maskstr[i] == '1'])
            fasta_outstr = ">%s\n%s\n" % (k, aligned_seq)
            outfwriter.write(fasta_outstr)
        outfwriter.close()
        return True

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
            self.logger.info("Check file {} ({}): OK".format(file_name, file_path))
            return True
        else:
            self.logger.warning("Check file {} ({}): missing".format(file_name, file_path))
            return False

    def checkfolder(self, folder_path, folder_name):
        """Check that a folder exists, output the result to the logger.

        Returns
        -------
        bool
            True if the folder exists, False otherwise.
        """
        if os.path.isdir(folder_path) and len(os.listdir(folder_path)) > 0:
            self.logger.info("Check folder {} ({}): OK".format(folder_name, folder_path))
            return True
        else:
            self.logger.warning("Check folder {} ({}): missing".format(folder_name, folder_path))
            return False

    def check_install(self):
        """Check that all reference files exist.

        Returns
        -------
        bool
            True if the installation is complete, False otherwise.
        """
        ok = True
        ok = ok and self.checkfile(Config.TAXONOMY_FILE, 'Taxonomy')
        ok = ok and self.checkfile(Config.CONCAT_BAC120, 'concat_bac120')
        ok = ok and self.checkfile(Config.CONCAT_AR122, 'concat_ar122')
        ok = ok and self.checkfile(os.path.join(Config.MASK_DIR, Config.MASK_BAC120), 'mask_bac120')
        ok = ok and self.checkfile(os.path.join(Config.MASK_DIR, Config.MASK_AR122), 'mask_ar122')
        ok = ok and self.checkfile(Config.TIGRFAM_HMMS, 'tirgfam_hmms')
        ok = ok and self.checkfile(os.path.join(Config.PFAM_HMM_DIR, 'Pfam-A.hmm'), 'pfam_hmms')

        ok = ok and self.checkfolder(Config.FASTANI_GENOMES, 'fastani_genomes')
        ok = ok and self.checkfolder(os.path.join(Config.PPLACER_DIR, Config.PPLACER_BAC120_REF_PKG), 'pplacer_bac120')
        ok = ok and self.checkfolder(os.path.join(Config.PPLACER_DIR, Config.PPLACER_AR122_REF_PKG), 'pplacer_ar122')

        if not ok:
            self.logger.error('One or more reference files are malformed.')
            raise ReferenceFileMalformed
