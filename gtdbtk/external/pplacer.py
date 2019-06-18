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
import subprocess

from gtdbtk.exceptions import PplacerException, TogException


class Pplacer(object):
    """ Phylogenetic placement of genomes into a tree (http://matsen.fredhutch.org/pplacer/). """

    def __init__(self):
        """ Instantiate the class. """
        self.logger = logging.getLogger('timestamp')

    def run(self, cpus, model, ref_pkg, json_out, msa_file, pplacer_out, mmap_file=None):
        """ Place genomes into the tree.

        Args:
            cpus (int): The number of threads to use.
            model (str): The model to use. Protein: LG, WAG, JTT. Nucleotides: GTR.
            ref_pkg (str): The path to the reference package.
            json_out (str): The path to write the json output to.
            msa_file (str): The path to the input MSA file.
            pplacer_out (str): Where to write the pplacer output file.
            mmap_file (str, optional): The path to write a scratch file to.

        Raises:
            PplacerException: if a non-zero exit code, or if the json output file isn't generated.

        """

        args = ['pplacer', '-m', model, '-j', str(cpus), '-c', ref_pkg, '-o', json_out, msa_file]

        if mmap_file:
            args.append('--mmap-file')
            args.append(mmap_file)

        with open(pplacer_out, 'w') as f_out:
            proc = subprocess.Popen(args, stdout=f_out, stderr=subprocess.PIPE)
            proc_out, proc_err = proc.communicate()

        if proc.returncode != 0:
            self.logger.error('An error was encountered while running pplacer.')
            raise PplacerException(proc_err)

        if not os.path.isfile(json_out):
            self.logger.error('pplacer returned a zero exit code but no output file was generated.')
            raise PplacerException

    def tog(self, pplacer_json_out, tree_file):
        """ Convert the pplacer json output into a newick tree.

        Args:
            pplacer_json_out (str): The path to the output of pplacer.
            tree_file (str): The path to output the newick file to.

        Raises:
            TogException: If a non-zero exit code is returned, or the tree file isn't output.
        """

        args = ['guppy', 'tog', '-o', tree_file, pplacer_json_out]
        proc = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        proc_out, proc_err = proc.communicate()

        if proc.returncode != 0:
            self.logger.error('An error was encountered while running tog.')
            raise TogException(proc_err)

        if not os.path.isfile(pplacer_json_out):
            self.logger.error('tog returned a zero exit code but no output file was generated.')
            raise TogException
