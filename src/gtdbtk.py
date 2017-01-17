#!/usr/bin/env python
###############################################################################
#                                                                             #
#                                                                             #
#                                                                             #
#    Entry point. See gtdbtk/main.py for internals                            #
#                                                                             #
#                                                                             #
#                                                                             #
###############################################################################
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

__author__ = "Pierre Chaumeil"
__copyright__ = "Copyright 2016-2017"
__credits__ = ["Pierre Chaumeil"]
__license__ = "GPL3"
__version__ = "0.0.1"
__maintainer__ = "Pierre Chaumeil"
__email__ = "uqpchaum@uq.edu.au"
__status__ = "Development"

###############################################################################

import argparse
import sys
from gtdbtk import gtdbtk

###############################################################################
###############################################################################
###############################################################################
###############################################################################


def printHelp():
    print '''\
                             ...::: gtdb toolkit :::...
   -------------------------------------------------------------------------
                                  version: %s
   -------------------------------------------------------------------------
    Typical workflow:
    gtdbtk tree        -> generate a Tree based on a multi sequence alignement.
    gtdbtk predict     -> predict if the user genomes are Archaeal or Bacterial genomes.
    USE: gtdbtk OPTION -h to see detailed options
    ''' % __version__

#    groopm import       -> Import data from csv


class CustomHelpFormatter(argparse.HelpFormatter):

    def _split_lines(self, text, width):
        return text.splitlines()

    def _get_help_string(self, action):
        h = action.help
        if '%(default)' not in action.help:
            if action.default != '' and \
               action.default != [] and \
               action.default is not None \
               and action.default != False:
                if action.default is not argparse.SUPPRESS:
                    defaulting_nargs = [argparse.OPTIONAL,
                                        argparse.ZERO_OR_MORE]

                    if action.option_strings or action.nargs in defaulting_nargs:

                        if '\n' in h:
                            lines = h.splitlines()
                            lines[0] += ' (default: %(default)s)'
                            h = '\n'.join(lines)
                        else:
                            h += ' (default: %(default)s)'
        return h

    def _fill_text(self, text, width, indent):
        return ''.join([indent + line for line in text.splitlines(True)])

if __name__ == '__main__':

    #-------------------------------------------------
    # intialise the options parser
    parser = argparse.ArgumentParser(prog='gtdb', add_help=False)
    parser.add_argument('-t', '--threads', type=int, default=1, help="Maximum number of threads/cpus to use.")
    parser.add_argument('-f', '--force', action="store_true", default=False, help="overwrite existing DB file without prompting")

    subparsers = parser.add_subparsers(help="--", dest='subparser_name')

    ##################################################
    # Typical workflow
    ##################################################

    #-------------------------------------------------
    # parse raw data and save
    tree_parser = subparsers.add_parser('tree',
                                        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                        help='generate a phylogenetic tree',)
    tree_parser.add_argument('--batchfile', help="Path to the User genomes")

    mutual_genome_add = tree_parser.add_argument_group('mutually exclusive required arguments')
    mutex_group = mutual_genome_add.add_mutually_exclusive_group(required=True)
    mutex_group.add_argument('--bacteria', action='store_true', dest='bac_domain')
    mutex_group.add_argument('--archaea', action='store_true', dest='arc_domain')

    tree_parser.add_argument('--no_tree', dest='no_tree', action="store_true",
                             help="Output tree data, but do not infer a tree.")
    tree_parser.add_argument('--min_perc_aa', type=float, default=50,
                             help='Filter genomes with an insufficient percentage of AA in the MSA.')
    tree_parser.add_argument('--consensus', type=float, default=25,
                             help='minimum percentage of the same amino acid required to retain column.')
    tree_parser.add_argument('--min_perc_taxa', type=float, default=50,
                             help='minimum percentage of taxa required required to retain column.')
    tree_parser.add_argument('--filter_taxa',
                             help='Filter genomes appearing on the output tree based on their toxonomic ranks(comma delimited).')
    tree_parser.add_argument('--prefix', required=False, default='gtdbtk',
                             help='Desired prefix for output files.')

    tree_parser.add_argument('--output', dest='out_dir', required=True,
                             help='Directory to output files.')

    #-------------------------------------------------
    # predict if the user genomes are Archaea or Bacteria
    predict_parser = subparsers.add_parser('predict',
                                           formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                           help='Predict if the user genomes are Archaeal or Bacterial genomes.')
    predict_parser.add_argument('--batchfile', help="Path to the User genomes")
    predict_parser.add_argument('--output_file', dest='outfile', help="Path to the output file")

    ##################################################
    # System
    ##################################################

    #-------------------------------------------------
    # get and check options
    args = None
    if(len(sys.argv) == 1):
        printHelp()
        sys.exit(0)
    elif(sys.argv[1] == '-v' or
         sys.argv[1] == '--v' or
         sys.argv[1] == '-version' or
         sys.argv[1] == '--version'):
        print "gtdbtk: version %s %s %s" % (__version__,
                                            __copyright__,
                                            __author__)
        sys.exit(0)
    elif(sys.argv[1] == '-h' or
         sys.argv[1] == '--h' or
         sys.argv[1] == '-help' or
         sys.argv[1] == '--help'):
        printHelp()
        sys.exit(0)
    else:
        args = parser.parse_args()

    #-------------------------------------------------
    # do what we came here to do
    try:
        GT_parser = gtdbtk.GtdbTKOptionsParser(__version__)
        if(False):
            import cProfile
            cProfile.run('GT_parser.parseOptions(args)', 'prof')
        else:
            GT_parser.parseOptions(args)
    except:
        print "Unexpected error:", sys.exc_info()[0]
        raise

###############################################################################
###############################################################################
###############################################################################
###############################################################################
