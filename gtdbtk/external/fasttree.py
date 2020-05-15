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

from gtdbtk.biolib_lite.common import make_sure_path_exists
from gtdbtk.biolib_lite.execute import check_dependencies
from gtdbtk.exceptions import FastTreeException


class FastTree(object):
    """Python wrapper for FastTree: http://www.microbesonline.org/fasttree/"""

    def __init__(self):
        """Instantiate the class."""
        self.logger = logging.getLogger('timestamp')
        self.version = self.get_version()

    def get_version(self):
        """ Get FastTree version. """
        try:
            env = os.environ.copy()
            args = ['FastTree']
            proc = subprocess.Popen(args, stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE, env=env, encoding='utf-8')
            output, error = proc.communicate()
            return error.split('\n')[0].split(' ')[4]
        except:
            return "(version unavailable)"

    def run(self, output_tree, tree_log, fasttree_log, prot_model, no_support, gamma, msa_file, cpus=1):
        """Run FastTree.

        Parameters
        ----------
        output_tree : str
            The path where the resulting tree should be written to.
        tree_log : str
            The path where the FastTree stats should be written to.
        fasttree_log : str
            The path where the FastTree log should be written to.
        prot_model : str
            Either 'JTT', 'WAG', or 'LG'.
        no_support : bool
            True if no support should be used, False otherwise.
        gamma : bool
            True if Gamma20 should be used, False otherwise.
        msa_file : str
            The path to the input MSA.
        cpus : int
            The maximum number of CPUs for FastTree to use.

        Raises
        ------
        FastTreeException
            If an error is encountered while running FastTree.

        """
        env = os.environ.copy()
        if cpus > 1:
            cmd = 'FastTreeMP'
            env['OMP_NUM_THREADS'] = str(cpus)
        else:
            cmd = 'FastTree'
        check_dependencies([cmd])

        make_sure_path_exists(os.path.dirname(output_tree))
        make_sure_path_exists(os.path.dirname(tree_log))
        make_sure_path_exists(os.path.dirname(fasttree_log))

        # Setup arguments
        args = [cmd]
        model_out = [prot_model]
        if prot_model == 'WAG':
            args.append('-wag')
        elif prot_model == 'LG':
            args.append('-lg')

        if gamma:
            args.append('-gamma')
            model_out.append('+G')

        if no_support:
            args.append('-nosupport')
        else:
            model_out.append('SH support values')

        args.append('-log')
        args.append(tree_log)
        args.append(msa_file)

        self.logger.info('Inferring FastTree ({}) using a maximum of {} CPUs.'.format(
            ', '.join(model_out), cpus))

        with open(output_tree, 'w') as f_out_tree:
            with open(fasttree_log, 'w') as f_out_err:
                proc = subprocess.Popen(
                    args, stdout=f_out_tree, stderr=f_out_err, env=env)
                proc.communicate()

        # Validate results
        if proc.returncode != 0:
            self.logger.error(
                'An error was encountered while running FastTree.')
            raise FastTreeException('FastTree returned a non-zero exit code.')
        if not os.path.isfile(output_tree):
            self.logger.error(
                'An error was encountered while running FastTree.')
            raise FastTreeException(
                'Tree output file is missing: {}'.format(output_tree))
        elif os.path.getsize(output_tree) < 1:
            self.logger.error(
                'An error was encountered while running FastTree.')
            raise FastTreeException(
                'Tree output file is empty: {}'.format(output_tree))
