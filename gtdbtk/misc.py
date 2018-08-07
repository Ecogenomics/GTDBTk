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

import os
import sys
import logging

import config.config as Config

from biolib.seq_io import read_fasta


class Misc():
    
    def __init__(self):
        """Initialize."""
        self.logger = logging.getLogger('timestamp') 

    def trim_msa(self,untrimmed_msa,mask_type,maskid,output_file):
        if maskid == 'bac' and mask_type == 'reference':
            mask = os.path.join(Config.MASK_DIR,Config.MASK_BAC120)
        elif maskid == 'arc' and mask_type == 'reference':
            mask = os.path.join(Config.MASK_DIR,Config.MASK_AR122)
        elif mask_type == 'file':
            mask = maskid
        with open(mask, 'r') as f:
            maskstr = f.readline()
            
        outfwriter = open(output_file,'w')
        dict_genomes = read_fasta(untrimmed_msa,False)
        
        for k,v in dict_genomes.iteritems():
            aligned_seq = ''.join([v[i] for i in xrange(0, len(maskstr)) if maskstr[i]=='1'])
            fasta_outstr = ">%s\n%s\n" % (k, aligned_seq)
            outfwriter.write(fasta_outstr)
        outfwriter.close()
        return True
    
    def checkfile(self,file_path,file_name):
        self.logger.warning("Check file {}: {}".format(file_name,file_path))
        if os.path.exists(file_path) and os.path.getsize(file_path) > 0:
            self.logger.warning("{} file..........OK".format(file_name))
        else:
            self.logger.warning("{} file..........missing".format(file_name))
            raise Exception("GTDB-Tk installation is incomplete.")
        
    def checkfolder(self,folder_path,folder_name):
        self.logger.warning("Check folder {}: {}".format(folder_name,folder_path))
        if os.path.isdir(folder_path) and len(os.listdir(folder_path)) > 0:
            self.logger.warning("{} dir..........OK".format(folder_name))
        else:
            self.logger.warning("{} dir..........missing".format(folder_name))
            raise Exception("GTDB-Tk installation is incomplete.")
    
    def check_install(self):
        try:
            self.checkfile(Config.TAXONOMY_FILE,'Taxonomy')
            self.checkfile(Config.CONCAT_BAC120,'concat_bac120')
            self.checkfile(Config.CONCAT_AR122,'concat_ar122')
            self.checkfile(os.path.join(Config.MASK_DIR,Config.MASK_BAC120),'mask_bac120')
            self.checkfile(os.path.join(Config.MASK_DIR,Config.MASK_AR122),'mask_ar122')
            self.checkfile(Config.TIGRFAM_HMMS,'tirgfam_hmms')
            pfam_test_file = os.path.join(Config.PFAM_HMM_DIR,'Pfam-A.hmm')
            self.checkfile(pfam_test_file,'pfam_hmms')
            
            self.checkfolder(Config.FASTANI_GENOMES,'fastani_genomes')
            self.checkfolder(os.path.join(Config.PPLACER_DIR,Config.PPLACER_BAC120_REF_PKG),'pplacer_bac120')
            self.checkfolder(os.path.join(Config.PPLACER_DIR,Config.PPLACER_AR122_REF_PKG),'pplacer_ar122')
        except Exception as e:
            raise
