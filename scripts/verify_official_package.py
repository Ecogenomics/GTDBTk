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

__prog_name__ = 'verify_official_package.py'
__prog_desc__ = 'Verify if all the files are properly formatted for an official release and archive the folder.'

__author__ = 'Pierre Chaumeil'
__copyright__ = 'Copyright 2019'
__credits__ = ['Pierre Chaumeil']
__license__ = 'GPL3'
__version__ = '0.0.2'
__maintainer__ = 'Pierre Chaumeil'
__email__ = 'p.chaumeil@uq.edu.au'
__status__ = 'Development'

import argparse
import copy
import glob
import os
import sys
import tarfile

import dendropy

from gtdbtk.biolib_lite.seq_io import read_fasta


class PackageChecker(object):
    def __init__(self, pack_dir):
        """Initialization."""
        self.pack_dir = pack_dir
        self.list_dirsinpackage = [
            'fastani', 'markers', 'masks', 'metadata', 'msa', 'pplacer', 'radii', 'taxonomy']

    def run(self, outf):

        # Check if all directories are here
        actual_dirs = os.listdir(self.pack_dir)
        if len(actual_dirs) != len(self.list_dirsinpackage):
            print('ERROR:')
        if len(set(actual_dirs) & set(self.list_dirsinpackage)) != len(self.list_dirsinpackage):
            print('ERROR:')

        with open(os.path.join(self.pack_dir, 'metadata', 'metadata.txt')) as metafile:
            for line in metafile:
                if line.startswith('VERSION_DATA'):
                    version = line.strip().split('=')[1]

        # List genomes in fastani folder
        list_genomes = [os.path.basename(x) for x in glob.glob(os.path.join(
            self.pack_dir, 'fastani', 'database/*.gz'))]
        list_genomes = [x.replace('_genomic.fna.gz', '').replace('GCA_', 'GB_GCA_').replace('GCF_', 'RS_GCF_') for x in
                        list_genomes]

        # Archaeal genome MSA is untrimmed
        ar_msa_file = glob.glob(os.path.join(
            self.pack_dir, 'msa/*ar122.faa'))[0]
        ar_msa = read_fasta(ar_msa_file)
        first_seq = ar_msa.get(list(ar_msa.keys())[0])
        if len(first_seq) != 32675:
            print('ERROR: len(first_seq) != 32675')

        # Bacterial genome MSA is untrimmed
        bac_msa_file = glob.glob(os.path.join(
            self.pack_dir, 'msa/*bac120.faa'))[0]
        bac_msa = read_fasta(bac_msa_file)
        first_seq = bac_msa.get(list(bac_msa.keys())[0])
        if len(first_seq) != 41155:
            print('ERROR: len(first_seq) != 41155')

        # Bacterial MASK is same length as the untrimmed bacterial genomes
        bac_mask_file = glob.glob(os.path.join(
            self.pack_dir, 'masks/*bac120.mask'))[0]
        bac_mask = ''
        with open(bac_mask_file) as bmf:
            bac_mask = bmf.readline()
        if len(bac_mask) != 41155:
            print('ERROR: len(bac_mask) != 41155')

        # Archaeal MASK is same length as the untrimmed archaeal genomes
        ar_mask_file = glob.glob(os.path.join(
            self.pack_dir, 'masks/*ar122.mask'))[0]
        ar_mask = ''
        with open(ar_mask_file) as amf:
            ar_mask = amf.readline()
        if len(ar_mask) != 32675:
            print('ERROR: len(ar_mask) != 32675')

        # Archaeal Pplacer MSA should have the same number of genomes as the
        # Archaeal untrimmed MSA
        ar_pplacer_msa_file = glob.glob(os.path.join(
            self.pack_dir, 'pplacer', 'gtdb_' + version + '_ar122.refpkg', 'ar122_msa_r95.faa'))[0]
        ar_pplacer_msa = read_fasta(ar_pplacer_msa_file)
        if len(ar_pplacer_msa) != len(ar_msa):
            print('ERROR: len(ar_pplacer_msa) != len(ar_msa)')
            print('len(ar_pplacer_msa): {}'.format(len(ar_pplacer_msa)))
            print('len(ar_msa): {}'.format(len(ar_msa)))
            print('difference genomes: {}'.format(list(set(ar_msa.keys()).difference(set(ar_pplacer_msa.keys())))))
        first_seq = ar_pplacer_msa.get(list(ar_pplacer_msa.keys())[0])
        # Archaeal Pplacer MSA should have the same length as the Archaeal mask
        if len(first_seq) != len([a for a in ar_mask if a == '1']):
            print('ERROR: len(first_seq) != len([a for a in ar_mask if a ==1])')
            print('len(first_seq): {}'.format(len(first_seq)))
            print('len([a for a in ar_mask if a ==1]): {}'.format(len([a for a in ar_mask if a == '1'])))

        # Bacterial Pplacer MSA should have the same number of genomes as the
        # Bacterial untrimmed MSA
        bac_pplacer_msa_file = os.path.join(
            self.pack_dir, 'pplacer', 'gtdb_' + version + '_bac120.refpkg', 'bac120_msa_r95.faa')
        bac_pplacer_msa = read_fasta(bac_pplacer_msa_file)
        if len(bac_pplacer_msa) != len(bac_msa):
            print('ERROR: len(bac_pplacer_msa) != len(bac_msa)')
            print('len(bac_pplacer_msa): {}'.format(len(bac_pplacer_msa)))
            print('len(bac_msa): {}'.format(len(bac_msa)))
            print('difference genomes: {}'.format(list(set(bac_msa.keys()).difference(set(bac_pplacer_msa.keys())))))
        first_seq = bac_pplacer_msa.get(list(bac_pplacer_msa.keys())[0])
        # Bacterial Pplacer MSA should have the same length as the Bacterial
        # mask
        if len(first_seq) != len([a for a in bac_mask if a == '1']):
            print('ERROR: len(first_seq) != len([a for a in bac_mask if a ==1])')
            print('len(first_seq): {}'.format(len(first_seq)))
            print('len([a for a in bac_mask if a ==1]): {}'.format(len([a for a in bac_mask if a == '1'])))

        # Archaeal Tree should have the same number of leaves than nomber of
        # genomes in the MSA
        arc_tree = dendropy.Tree.get_from_path(os.path.join(
            self.pack_dir, 'pplacer', 'gtdb_' + version + '_ar122.refpkg', 'ar122_' + version + '_unroot.pplacer.tree'),
            schema='newick',
            rooting='force-rooted',
            preserve_underscores=True)
        list_leaves = arc_tree.leaf_nodes()
        if len(list_leaves) != len(ar_pplacer_msa):
            print('ERROR: len(list_leaves) != len(ar_pplacer_msa)')
            print('len(list_leaves): {}'.format(len(list_leaves)))
            print('len(ar_pplacer_msa): {}'.format(len(ar_pplacer_msa)))

        # Bacterial Tree should have the same number of leaves than nomber of
        # genomes in the MSA
        bac_tree = dendropy.Tree.get_from_path(os.path.join(
            self.pack_dir, 'pplacer', 'gtdb_' + version + '_bac120.refpkg',
                                      'bac120_' + version + '_unroot.pplacer.tree'),
            schema='newick',
            rooting='force-rooted',
            preserve_underscores=True)
        list_leaves = bac_tree.leaf_nodes()
        if len(list_leaves) != len(bac_pplacer_msa):
            print('ERROR: len(list_leaves) != len(bac_pplacer_msa)')
            print('len(list_leaves): {}'.format(len(list_leaves)))
            print('len(bac_pplacer_msa): {}'.format(len(bac_pplacer_msa)))

        # Taxonomy file should have as many genomes as bac120 and ar122 MSA
        # combined
        tax_file = os.path.join(
            self.pack_dir, 'taxonomy', 'gtdb_taxonomy.tsv')
        tax_dict = {}
        with open(tax_file) as tf:
            for line in tf:
                infos = line.strip().split('\t')
                tax_dict[infos[0]] = infos[1]
        if len(tax_dict) != (len(ar_msa) + len(bac_msa)):
            print('ERROR: len(tax_dict) != (len(ar_msa) + len(bac_msa))')
            print('len(tax_dict): {}'.format(len(tax_dict)))
            print('len(ar_msa) + len(bac_msa): {}'.format(len(ar_msa) + len(bac_msa)))

        # Radii file should have as many genomes as bac120 and ar122 MSA
        # combined
        radii_file = os.path.join(
            self.pack_dir, 'radii', 'gtdb_radii.tsv')
        radii_dict = {}
        with open(radii_file) as rf:
            for line in rf:
                infos = line.strip().split('\t')
                radii_dict[infos[1]] = infos[2]
        if len(radii_dict) != (len(ar_msa) + len(bac_msa)):
            print('ERROR: len(radii_dict) != (len(ar_msa) + len(bac_msa))')
            print('len(radii_dict): {}'.format(len(radii_dict)))
            print('len(ar_msa) + len(bac_msa): {}'.format(len(ar_msa) + len(bac_msa)))
        if len(set(radii_dict.keys()).symmetric_difference(set(tax_dict.keys()))) != 0:
            print('ERROR: len(set(radii_dict.keys()).symmetric_difference(tax_dict.keys()))')
            print('set(radii_dict.keys()).symmetric_difference(tax_dict.keys()): {}'.format(
                set(radii_dict.keys()).symmetric_difference(set(tax_dict.keys()))))

        if len(list_genomes) != len(radii_dict):
            print('ERROR: len(list_genomes) != len(radii_dict)')
            print('Missing genomes {}'.format(set(list_genomes) ^ set(radii_dict.keys())))
            print('len(list_genomes): {}'.format(len(list_genomes)))
            print('len(radii_dict): {}'.format(len(radii_dict)))

        print('\n\nVERSION: {}'.format(version))
        print('Length trimmed bac120 MSA: {}'.format(len(bac_pplacer_msa.get(list(bac_pplacer_msa.keys())[0]))))
        print('Length trimmed ar122 MSA: {}'.format(len(ar_pplacer_msa.get(list(ar_pplacer_msa.keys())[0]))))
        print('')
        print('Number of genomes in fastani/database: {}'.format(len(list_genomes)))
        print('Number of genomes in radii file: {}'.format(len(radii_dict)))
        print('Number of genomes in taxonomy file: {}'.format(len(tax_dict)))

        print('Would you like to archive the folder? ')
        # raw_input returns the empty string for "enter"

        yes = {'yes', 'y', 'yep', ''}
        no = {'no', 'n'}

        final_choice = False
        choice = input().lower()
        if choice in yes:
            with tarfile.open(outf, "w:gz") as tar:
                packdir = copy.copy(self.pack_dir)
                if packdir.endswith('/'):
                    packdir = packdir[:-1]
                tar.add(self.pack_dir, arcname=os.path.basename(packdir))
        elif choice in no:
            return False
        else:
            sys.stdout.write("Please respond with 'yes' or 'no'")


if __name__ == '__main__':
    print(__prog_name__ + ' v' + __version__ + ': ' + __prog_desc__)
    print('  by ' + __author__ + ' (' + __email__ + ')' + '\n')

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--package_directory',
                        help='Directory where the unarchived package is located.')
    parser.add_argument(
        '--output_file', help='File name similar to gtdbtk.r86_v3_data.tar.gz .')
    args = parser.parse_args()

    try:
        package_checker = PackageChecker(args.package_directory)
        package_checker.run(args.output_file)
    except SystemExit:
        print("\nControlled exit resulting from an unrecoverable error or warning.")
    except:
        print("\nUnexpected error:", sys.exc_info()[0])
        raise
