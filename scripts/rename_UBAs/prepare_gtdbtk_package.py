#!/usr/bin/env python

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

__prog_name__ = 'prepare_gtdbtk_package.py'
__prog_desc__ = 'Rename user genomes to their GCAs or UBAs ids.'

__author__ = 'Pierre Chaumeil'
__copyright__ = 'Copyright 2018'
__credits__ = ['Pierre Chaumeil']
__license__ = 'GPL3'
__version__ = '0.0.1'
__maintainer__ = 'Pierre Chaumeil'
__email__ = 'uqpchaum@uq.edu.au'
__status__ = 'Development'

import argparse
import glob
import logging
import os
import sys
from shutil import copyfile

import dendropy

from gtdbtk.biolib_lite.custom_help_formatter import CustomHelpFormatter
from gtdbtk.biolib_lite.logger import logger_setup
from gtdbtk.biolib_lite.seq_io import read_fasta


class OptionsParser(object):

    def __init__(self):
        """Initialization"""

        self.logger = logging.getLogger('timestamp')

    def parse_options(self, options):
        """Parse user options and call the correct pipeline(s)"""

        if options.subparser_name == 'convert_skani_genomes':
            self.convert_skani_genomes(options.dirin, options.dirout)
        elif options.subparser_name == 'unroot_tree':
            self.unroot(options.input_tree, options.output_tree)
        else:
            self.logger.error('Unknown command: ' +
                              options.subparser_name + '\n')
            sys.exit()

        return 0

    def convert_skani_genomes(self, dirin, dirout):

        # get mapping from published UBA genomes to NCBI accessions
        __location__ = os.path.realpath(os.path.join(
            os.getcwd(), os.path.dirname(__file__)))

        uba_acc = {}
        with open(os.path.join(__location__, 'uba_ncbi_accessions.tsv')) as ub:
            for line in ub:
                line_split = line.strip('\n').split('\t')
                if line_split[2].startswith('G'):
                    uba_acc[line_split[0]] = {
                        "uba": line_split[1], "gca": '' + line_split[2]}
                else:
                    uba_acc[line_split[0]] = {"uba": line_split[1]}

        skanis = glob.glob(os.path.join(dirin, "*"))
        skani_dir = os.path.join(dirout, 'skani')
        if not os.path.exists(skani_dir):
            os.makedirs(skani_dir)
        for genome in skanis:
            filenamef = os.path.basename(genome)
            filenamef = filenamef.replace("_genomic.fna", "")
            if filenamef.startswith("U_"):
                print(filenamef)
                subdict = uba_acc.get(filenamef)
                if "gca" in subdict.keys():
                    copyfile(genome, os.path.join(
                        skani_dir, subdict.get("gca")[3:] + "_genomic.fna"))
                else:
                    copyfile(genome, os.path.join(
                        skani_dir, subdict.get("uba") + "_genomic.fna"))
            else:
                copyfile(genome, os.path.join(
                    skani_dir, filenamef + "_genomic.fna"))

    def unroot(self, rooted_tree, unrooted_tree):
        # get taxonomic classification of each user genome
        tree = dendropy.Tree.get_from_path(rooted_tree,
                                           schema='newick',
                                           rooting='force-rooted',
                                           preserve_underscores=True)
        tree.deroot()
        tree.write_to_path(unrooted_tree,
                           schema='newick',
                           suppress_rooting=True,
                           unquoted_underscores=True)
        return True

    def run(self, dirin, dirout, gtr, release):
        """ renaming genome files for skani"""

        # get list of genomes to retain (based on genome list 1014)
        genomes_to_retain = set()
        with open(gtr) as f:
            # f.readline()

            for line in f:
                line_split = line.strip().split('\t')
                genomes_to_retain.add(line_split[0])

        print('Genome to retain: %d' % len(genomes_to_retain))
        # get mapping from published UBA genomes to NCBI accessions
        __location__ = os.path.realpath(os.path.join(
            os.getcwd(), os.path.dirname(__file__)))

        uba_acc = {}
        with open(os.path.join(__location__, 'uba_ncbi_accessions.tsv')) as ub:
            for line in ub:
                line_split = line.strip().split('\t')
                if line_split[2] != "None":
                    uba_acc[line_split[0]] = {
                        "uba": line_split[1], "gca": 'GB_' + line_split[2]}
                else:
                    uba_acc[line_split[0]] = {"uba": line_split[1]}

        # renaming taxonomy:
        taxout = open(os.path.join(dirout, 'gtdb_taxonomy.tsv'), 'w')
        with open(os.path.join(dirin, 'gtdb_taxonomy.tsv')) as gt:
            for line in gt:
                info = line.strip().split("\t")
                if info[0] in genomes_to_retain:
                    if info[0].startswith("U_"):
                        subdict = uba_acc.get(info[0])
                        if "gca" in subdict.keys():
                            taxout.write("{0}\t{1}\n".format(
                                subdict.get("gca"), info[1]))
                        else:
                            taxout.write("{0}\t{1}\n".format(
                                subdict.get("uba"), info[1]))
                    else:
                        taxout.write(line)
        taxout.close()

        # renaming genome files for skani
        skanis = glob.glob(os.path.join(dirin, 'skani', "*"))
        skani_dir = os.path.join(dirout, 'skani')
        if not os.path.exists(skani_dir):
            os.makedirs(skani_dir)
        for genome in skanis:
            filenamef = os.path.basename(genome)
            filenamef = filenamef.replace("_genomic.fna", "")
            if filenamef.startswith("U_"):
                subdict = uba_acc.get(filenamef)
                if filenamef == "U_74684":
                    print(subdict)
                    print(genome)
                    print(os.path.join(skani_dir, subdict.get("gca")[3:] + "_genomic.fna"))
                if "gca" in subdict.keys():
                    copyfile(genome, os.path.join(
                        skani_dir, subdict.get("gca")[3:] + "_genomic.fna"))
                else:
                    copyfile(genome, os.path.join(
                        skani_dir, subdict.get("uba") + "_genomic.fna"))
            else:
                copyfile(genome, os.path.join(
                    skani_dir, filenamef + "_genomic.fna"))

        for dom in ['bac120', 'ar53']:
            # MSA renaming
            msadir = os.path.join(dirout, dom, 'msa')
            if not os.path.exists(msadir):
                os.makedirs(msadir)
            msa_dict = read_fasta(os.path.join(
                dirin, dom, 'gtdb_concatenated.faa'))
            seqout = open(os.path.join(msadir, 'gtdb_r' +
                                       release + '_' + dom + '.faa'), 'w')
            for gid, seq in msa_dict.items():
                if gid in genomes_to_retain:
                    if gid.startswith("U_"):
                        subdict = uba_acc.get(gid)
                        if "gca" in subdict.keys():
                            seqout.write(">{0}\n{1}\n".format(
                                subdict.get("gca"), seq))
                        else:
                            seqout.write(">{0}\n{1}\n".format(
                                subdict.get("uba"), seq))
                    else:
                        seqout.write(">{0}\n{1}\n".format(gid, seq))
            seqout.close()

            # PPLACER renaming
            pplacerdir = os.path.join(dirout, dom, 'pplacer')
            if not os.path.exists(pplacerdir):
                os.makedirs(pplacerdir)

            trees = glob.glob(os.path.join(dirin, dom, 'pplacer', "*.tree"))
            if len(trees) != 1:
                print("Error")
                sys.exit()
            else:
                treef = trees[0]
            fastas = glob.glob(os.path.join(dirin, dom, 'pplacer', "*.fa"))
            if len(fastas) != 1:
                print("Error")
                sys.exit()
            else:
                seqfile = fastas[0]
            logs = glob.glob(os.path.join(dirin, dom, 'pplacer', "*.log"))
            if len(logs) != 1:
                print("Error")
                sys.exit()
            else:
                logfile = logs[0]

            # produce corrected tree
            tree = dendropy.Tree.get_from_path(os.path.join(treef),
                                               schema='newick',
                                               rooting='force-rooted',
                                               preserve_underscores=True)
            for n in tree.leaf_node_iter():
                if n.taxon.label.startswith("U_"):
                    subdict = uba_acc.get(n.taxon.label)
                    if "gca" in subdict.keys():
                        n.taxon.label = subdict.get("gca")
                    else:
                        n.taxon.label = subdict.get("uba")
            tree.write_to_path(os.path.join(dirout, dom, 'pplacer', dom + "_r" + release + ".tree"),
                               schema='newick',
                               suppress_rooting=True,
                               unquoted_underscores=True)

            trimmed_seqout = open(os.path.join(
                dirout, dom, 'pplacer', 'trimmed_msa_' + dom + '.faa'), 'w')
            trimmed_fasta = read_fasta(seqfile)
            for gid, seq in trimmed_fasta.items():
                if gid in genomes_to_retain:
                    if gid.startswith("U_"):
                        subdict = uba_acc.get(gid)
                        if "gca" in subdict.keys():
                            trimmed_seqout.write(
                                ">{0}\n{1}\n".format(subdict.get("gca"), seq))
                        else:
                            trimmed_seqout.write(
                                ">{0}\n{1}\n".format(subdict.get("uba"), seq))
                    else:
                        trimmed_seqout.write(">{0}\n{1}\n".format(gid, seq))
            trimmed_seqout.close()

            logoutf = open(os.path.join(dirout, dom, 'pplacer',
                                        'fitting_' + dom + '.log'), 'w')
            with open(logfile) as logfin:
                for line in logfin:
                    for k, subdict in uba_acc.items():
                        if "gca" in subdict.keys():
                            line = line.replace(
                                k + ":", subdict.get("gca") + ":")
                        else:
                            line = line.replace(
                                k + ":", subdict.get("uba") + ":")
                    logoutf.write(line)
            logoutf.close()


def print_help():
    """Help menu."""

    print('')
    print('                ...::: GTDB Tk Toolset v' + __version__ + ' :::...''')
    print('''\

    Moving to new release:
      blah
      blah
      blah

    External files for GTDB-Tk package:
      convert_skani_genomes  -> Rename genomic.fna files to UBAs or GCAs.
      unroot_tree              -> Unroot tree.

    ''')


if __name__ == "__main__":
    print(__prog_name__ + ' v' + __version__ + ': ' + __prog_desc__)
    print('  by ' + __author__ + ' (' + __email__ + ')' + '\n')

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    subparsers = parser.add_subparsers(help="--", dest='subparser_name')

    # generate convert genomic file for GTDB website
    skani_convert_parser = subparsers.add_parser('convert_skani_genomes',
                                                   formatter_class=CustomHelpFormatter,
                                                   description='rename Genomic.fna files to UBAs or GCAs.')
    skani_convert_parser.add_argument('--input', dest="dirin",
                                        required=True, help='path to taxonomy file.')
    skani_convert_parser.add_argument('--output', dest="dirout",
                                        required=True, help='path to processed output file')
    skani_convert_parser.add_argument(
        '--silent', help="suppress output", action='store_true')

    # generate convert genomic file for GTDB website
    unroot_tree_parser = subparsers.add_parser('unroot_tree',
                                               formatter_class=CustomHelpFormatter,
                                               description='Unroot tree.')
    unroot_tree_parser.add_argument('input_tree', help='')
    unroot_tree_parser.add_argument('output_tree', help='')
    unroot_tree_parser.add_argument(
        '--silent', help="suppress output", action='store_true')

    # get and check options
    args = None
    if len(sys.argv) == 1 or sys.argv[1] == '-h' or sys.argv == '--help':
        print_help()
        sys.exit(0)
    else:
        args = parser.parse_args()

    try:
        logger_setup(args.output_dir,
                     'gtdbtk_toolset.log',
                     'GTDB Tk converter',
                     version(),
                     args.silent)
    except:
        logger_setup(None,
                     'gtdbtk_toolset.log',
                     'GTDB Tk converter',
                     __version__,
                     args.silent)

    # do what we came here to do
    try:
        parser = OptionsParser()
        if False:
            # import pstats
            # p = pstats.Stats('prof')
            # p.sort_stats('cumulative').print_stats(10)
            # p.sort_stats('time').print_stats(10)
            import cProfile

            cProfile.run('parser.parse_options(args)', 'prof')
        elif False:
            import pdb

            pdb.run(parser.parse_options(args))
        else:
            parser.parse_options(args)
    except SystemExit:
        print("\n  Controlled exit resulting from an unrecoverable error or warning.")
    except:
        print("\nUnexpected error:", sys.exc_info()[0])
        raise
