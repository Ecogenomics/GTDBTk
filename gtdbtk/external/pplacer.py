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
import re
import subprocess
import sys

from gtdbtk.exceptions import PplacerException, TogException


class PplacerLogger(object):
    """Helper class for writing pplacer output."""

    def __init__(self, fh):
        """Initialise the class.

        Parameters
        ----------
        fh : BinaryIO
            The file to write to .
        """
        self.fh = fh
        self.init_dots = 0

    def _disp_progress(self, line):
        """Calculates the progress and writes it to stdout.

        Parameters
        ----------
        line : str
            The line passed from pplacer stdout.
        """
        if not line.startswith('working on '):
            sys.stdout.write('\rInitialising pplacer{}'.format('.' *
                                                               self.init_dots))
            sys.stdout.flush()
            self.init_dots = (self.init_dots + 1) % 4
        else:
            re_hits = re.search(r'\((\d+)\/(\d+)\)', line)
            current = int(re_hits.group(1))
            total = int(re_hits.group(2))
            sys.stdout.write('\r{}'.format(self._get_progress_str(current,
                                                                  total)))
            sys.stdout.flush()

    def _get_progress_str(self, current, total):
        """Determines the format of the genomes % string.

        Parameters
        ----------
        current : int
            The current number of genomes which have been placed.
        total : int
            The total number of genomes which are to be placed.

        Returns
        -------
        out : str
            A string formatted to show the progress of placement.
        """
        width = 50
        bar = str()
        prop = float(current) / total
        bar += '#' * int(prop * width)
        bar += '-' * (width - len(bar))
        return 'Placing genomes |{}| {}/{} ({:.2f}%)'.format(bar, current,
                                                             total, prop * 100)

    def read(self, line):
        """Reads a line and writes the progress to stdout and the file.

        Parameters
        ----------
        line : str
            A line returned from Prodigal stdout.
        """
        self.fh.write(line)
        line = line.strip()
        self._disp_progress(line)


class Pplacer(object):
    """Phylogenetic placement of genomes into a reference tree
    (http://matsen.fredhutch.org/pplacer/).
    """

    def __init__(self):
        """ Instantiate the class. """
        self.logger = logging.getLogger('timestamp')

    def run(self, cpus, model, ref_pkg, json_out, msa_file, pplacer_out,
            mmap_file=None):
        """Place genomes into a reference tree.

        Args:
            cpus (int): The number of threads to use.
            model (str): The model to use. PROT: LG, WAG, JTT. NT: GTR.
            ref_pkg (str): The path to the reference package.
            json_out (str): The path to write the json output to.
            msa_file (str): The path to the input MSA file.
            pplacer_out (str): Where to write the pplacer output file.
            mmap_file (str, optional): The path to write a scratch file to.

        Raises:
            PplacerException: if a non-zero exit code, or if the json output
                              file isn't generated.

        """

        args = ['pplacer', '-m', model, '-j', str(cpus), '-c', ref_pkg, '-o',
                json_out, msa_file]

        if mmap_file:
            args.append('--mmap-file')
            args.append(mmap_file)

        proc = subprocess.Popen(args, stdout=subprocess.PIPE)
        with open(pplacer_out, 'w') as fh:
            pplacer_logger = PplacerLogger(fh)
            while True:
                line = proc.stdout.readline()
                if not line:
                    sys.stdout.write('\n')
                    break
                pplacer_logger.read(line)
        proc.wait()

        if proc.returncode != 0:
            raise PplacerException('An error was encountered while '
                                   'running pplacer, check the log '
                                   'file: {}'.format(pplacer_out))

        if not os.path.isfile(json_out):
            self.logger.error('pplacer returned a zero exit code but no output '
                              'file was generated.')
            raise PplacerException

    def tog(self, pplacer_json_out, tree_file):
        """ Convert the pplacer json output into a newick tree.

        Args:
            pplacer_json_out (str): The path to the output of pplacer.
            tree_file (str): The path to output the newick file to.

        Raises:
            TogException: If a non-zero exit code is returned, or the tree file
                          isn't output.
        """

        args = ['guppy', 'tog', '-o', tree_file, pplacer_json_out]
        proc = subprocess.Popen(args, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
        proc_out, proc_err = proc.communicate()

        if proc.returncode != 0:
            self.logger.error('An error was encountered while running tog.')
            raise TogException(proc_err)

        if not os.path.isfile(pplacer_json_out):
            self.logger.error('tog returned a zero exit code but no output '
                              'file was generated.')
            raise TogException
