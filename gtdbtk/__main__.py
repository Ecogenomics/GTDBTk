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
import sys
import traceback

from gtdbtk import __author__, __copyright__, __version__
from gtdbtk.cli import get_main_parser
from gtdbtk.biolib_lite.exceptions import BioLibError
from gtdbtk.biolib_lite.logger import logger_setup
from gtdbtk.exceptions import *
from gtdbtk.main import OptionsParser
from gtdbtk.tools import get_gtdbtk_latest_version


def print_help():
    print('''\

              ...::: GTDB-Tk v%s :::...

  Workflows:
    classify_wf -> Classify genomes by placement in GTDB reference tree
                     (ani_screening -> identify -> align -> classify)
    de_novo_wf  -> Infer de novo tree and decorate with GTDB taxonomy
                     (identify -> align -> infer -> root -> decorate)

  Methods:
    identify -> Identify marker genes in genome
    align    -> Create multiple sequence alignment
    classify -> Determine taxonomic classification of genomes
    infer    -> Infer tree from multiple sequence alignment
    root     -> Root tree using an outgroup
    decorate -> Decorate tree with GTDB taxonomy

  Tools:
    infer_ranks     -> Establish taxonomic ranks of internal nodes using RED
    ani_rep         -> Calculates ANI to GTDB representative genomes
    trim_msa        -> Trim an untrimmed MSA file based on a mask
    export_msa      -> Export the untrimmed archaeal or bacterial MSA file
    remove_labels   -> Remove labels (bootstrap values, node labels) from an Newick tree
    convert_to_itol -> Convert a GTDB-Tk Newick tree to an iTOL tree
 

  Testing:
    test          -> Validate the classify_wf pipeline with 3 archaeal genomes 
    check_install -> Verify third party programs and GTDB reference package

  Use: gtdbtk <command> -h for command specific help
    ''' % __version__)


def main():
    # -------------------------------------------------
    # get and check options
    args = None
    if len(sys.argv) == 1:
        print_help()
        sys.exit(0)
    elif sys.argv[1] in {'-v', '--v', '-version', '--version'}:
        print(f"gtdbtk: version {__version__} {__copyright__} {__author__}")

        # Warn the user they are not using the latest version (if possible)
        latest_ver = get_gtdbtk_latest_version()
        if latest_ver and latest_ver != __version__:
            print(f'There is a newer version of GTDB-Tk available: v{latest_ver}')
        sys.exit(0)
    elif sys.argv[1] in {'-h', '--h', '-help', '--help'}:
        print_help()
        sys.exit(0)
    else:
        args = get_main_parser().parse_args()

    # setup logger
    logger_setup(args.out_dir if hasattr(args, 'out_dir') and args.out_dir else None,
                 "gtdbtk.log", "GTDB-Tk", __version__, False,
                 hasattr(args, 'debug') and args.debug)
    logger = logging.getLogger('timestamp')


    # -------------------------------------------------
    # do what we came here to do
    try:
        gt_parser = OptionsParser(__version__,args.out_dir if hasattr(args, 'out_dir') and args.out_dir else None)
        gt_parser.parse_options(args)
    except SystemExit:
        logger.error('Controlled exit resulting from early termination.')
        sys.exit(1)
    except KeyboardInterrupt:
        logger.error('Controlled exit resulting from interrupt signal.')
        sys.exit(1)
    except GTDBTkExit as e:
        if len(str(e)) > 0:
            logger.error('{}'.format(e))
        logger.error('Controlled exit resulting from an unrecoverable error or warning.')
        sys.exit(1)
    except (GTDBTkException, BioLibError) as e:
        msg = 'Controlled exit resulting from an unrecoverable error or warning.\n\n'
        msg += '=' * 80 + '\n'
        msg += 'EXCEPTION: {}\n'.format(type(e).__name__)
        msg += '  MESSAGE: {}\n'.format(e)
        msg += '_' * 80 + '\n\n'
        msg += traceback.format_exc()
        msg += '=' * 80
        logger.error(msg)
        sys.exit(1)
    except Exception as e:
        msg = 'Uncontrolled exit resulting from an unexpected error.\n\n'
        msg += '=' * 80 + '\n'
        msg += 'EXCEPTION: {}\n'.format(type(e).__name__)
        msg += '  MESSAGE: {}\n'.format(e)
        msg += '_' * 80 + '\n\n'
        msg += traceback.format_exc()
        msg += '=' * 80
        logger.error(msg)
        sys.exit(1)


if __name__ == '__main__':
    main()
