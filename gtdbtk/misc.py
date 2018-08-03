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

import config.config as Config
import os
import sys
from biolib.seq_io import read_fasta


class Misc():
    @staticmethod
    def trim_msa(untrimmed_msa,mask_type,maskid,output_file):
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
