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

import shutil

import dendropy

import gtdbtk.config.config as Config
from gtdbtk.biolib_lite.execute import check_dependencies
from gtdbtk.biolib_lite.logger import colour
from gtdbtk.biolib_lite.newick import parse_label
from gtdbtk.biolib_lite.seq_io import read_fasta
from gtdbtk.config.output import DIR_CLASSIFY_INTERMEDIATE, DIR_ALIGN_INTERMEDIATE, DIR_IDENTIFY_INTERMEDIATE
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
            mask = os.path.join(Config.MASK_DIR, Config.MASK_AR53)
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

    def remove_labels(self, input_file, output_file):
        """Remove labels from a Newick Tree.

        Parameters
        ----------
        input_file : str
            The path to the input Newick tree.
        output_file : str
            The path to the output Newick tree.
        """

        self.logger.info("Removing labels from tree {}".format(input_file))
        intree= dendropy.Tree.get_from_path(input_file,
                                           schema='newick',
                                           rooting='force-rooted',
                                           preserve_underscores=True)

        for node in intree.internal_nodes():
            node.label = None

        intree.write_to_path(output_file, schema='newick', suppress_rooting=True,unquoted_underscores=True)


    def convert_to_itol(self, input_file, output_file):
        """Remove labels from a Newick Tree.

        Parameters
        ----------
        input_file : str
            The path to the input Newick tree.
        output_file : str
            The path to the output Newick tree.
        """

        self.logger.info("Convert GTDB-Tk tree to iTOL format")
        intree= dendropy.Tree.get_from_path(input_file,
                                           schema='newick',
                                           rooting='force-rooted',
                                           preserve_underscores=True)

        for node in intree.internal_nodes():
            if node.label:
                bootstrap,label,_aux = parse_label(node.label)
                if label:
                    label = label.replace('; ',';').replace(';','|').replace("'","").lstrip('')
                node.label = label
                if node.edge.length:
                    node.edge.length = f'{node.edge.length}[{bootstrap}]'

        intree.write_to_path(output_file, schema='newick', suppress_rooting=True,unquoted_underscores=True)


    def remove_intermediate_files(self,output_dir,wf_name):
        """Remove intermediate files.

        Parameters
        ----------
        output_dir : str
            The path to the output directory.
        wf_name : str
            The name of the workflow to delete intermediate files.
        """
        self.logger.info('Removing intermediate files.')
        #Remove identify step intermediate files
        intermediate_identify = os.path.join(output_dir, DIR_IDENTIFY_INTERMEDIATE)
        if os.path.exists(intermediate_identify) and os.path.isdir(intermediate_identify):
            shutil.rmtree(intermediate_identify)
        #Remove align step intermediate files
        intermediate_align = os.path.join(output_dir, DIR_ALIGN_INTERMEDIATE)
        if os.path.exists(intermediate_align) and os.path.isdir(intermediate_align):
            shutil.rmtree(intermediate_align)
        if wf_name == 'classify_wf':
            #Remove classify step intermediate files
            intermediate_classify = os.path.join(output_dir, DIR_CLASSIFY_INTERMEDIATE)
            if os.path.exists(intermediate_classify) and os.path.isdir(intermediate_classify):
                shutil.rmtree(intermediate_classify)
        elif wf_name == 'de_novo_wf':
            #Remove classify step intermediate files
            intermediate_infer = os.path.join(output_dir, DIR_ALIGN_INTERMEDIATE)
            if os.path.exists(intermediate_infer) and os.path.isdir(intermediate_infer):
                shutil.rmtree(intermediate_infer)
        self.logger.info('Intermediate files removed.')

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
                    base_name, colour(f'HASH MISMATCH {user_hash}', ['bright'], fg='yellow')))
                ok = False
            else:
                self.logger.info("         |-- {:16} {}".format(
                    base_name, colour('OK', ['bright'], fg='green')))

        if not ok:
            raise GTDBTkExit('Unexpected files were seen, or the reference package is corrupt.')
