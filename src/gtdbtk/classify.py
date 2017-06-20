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
import shutil
import logging
import tempfile

from biolib.common import remove_extension
from biolib.seq_io import read_seq, read_fasta
from biolib.newick import parse_label
from biolib.external.execute import check_dependencies
from biolib.taxonomy import Taxonomy

from tools import genomes_to_process,add_ncbi_prefix,merge_two_dicts


import config.config as Config

import dendropy


class Classify(object):
    """Determine taxonomic classification of genomes by ML placement."""

    def __init__(self, cpus=1):
        """Initialize."""
        
        check_dependencies(['pplacer', 'guppy'])
        
        self.taxonomy_file = Config.TAXONOMY_FILE

        
        self.logger = logging.getLogger('timestamp') 
        self.cpus = cpus

    def place_genomes(self, 
                        user_msa_file, 
                        marker_set_id, 
                        out_dir, 
                        prefix):
        """Place genomes into reference tree using pplacer."""

        # get path to pplacer reference package   
        if marker_set_id == 'bac120':
            pplacer_ref_pkg = os.path.join(Config.PPLACER_DIR, Config.PPLACER_BAC120_REF_PKG)
        elif marker_set_id == 'ar122':
            pplacer_ref_pkg = os.path.join(Config.PPLACER_DIR, Config.PPLACER_AR122_REF_PKG)
        elif marker_set_id == 'rps23':
            pplacer_ref_pkg = os.path.join(Config.PPLACER_DIR, Config.PPLACER_RPS23_REF_PKG)
            
        # rename user MSA file for compatibility with pplacer
        if not user_msa_file.endswith('.fasta'):
            t = os.path.join(out_dir, prefix + '.user_msa.fasta')
            shutil.copyfile(user_msa_file, t)
            user_msa_file = t
              
        # run pplacer to place bins in reference genome tree
        num_genomes = sum([1 for _seq_id, _seq in read_seq(user_msa_file)])
        self.logger.info('Placing %d genomes into GTDB reference tree with pplacer (be patient).' % num_genomes)
        
        pplacer_out_dir = os.path.join(out_dir, 'classify', 'pplacer')
        if not os.path.exists(pplacer_out_dir):
            os.makedirs(pplacer_out_dir)
            
        pplacer_out = os.path.join(pplacer_out_dir, Config.PPLACER_OUT)
        pplacer_json_out = os.path.join(pplacer_out_dir, Config.PPLACER_JSON_OUT)
        cmd = 'pplacer -j %d -c %s -o %s %s > %s' % (self.cpus,
                                                     pplacer_ref_pkg,
                                                     pplacer_json_out,
                                                     user_msa_file,
                                                     pplacer_out)
        
        os.system(cmd)

        # extract tree
        tree_file = os.path.join(out_dir, 'classify', prefix + ".classify.tree")
        cmd = 'guppy tog -o %s %s' % (tree_file, pplacer_json_out)
        os.system(cmd)
        
        return tree_file
        
    def run(self, 
            user_msa_file,
            genome_dir, 
            batchfile, 
            marker_set_id, 
            out_dir, 
            prefix):
        """Classify genomes based on position in reference tree."""
        
        #=======================================================================
        # classify_tree = self.place_genomes(user_msa_file,
        #                                     marker_set_id,
        #                                     out_dir,
        #                                     prefix)
        #=======================================================================
        
        genomes = genomes_to_process(genome_dir, batchfile)
        classify_tree = '/srv/home/uqpchaum/playground/gtdbtk_test/classify_10UBAs/classify/gtdbtk.classify.tree'
        # get taxonomic classification of each user genome
        tree = dendropy.Tree.get_from_path(classify_tree, 
                                            schema='newick', 
                                            rooting='force-rooted', 
                                            preserve_underscores=True)
        
        gtdb_taxonomy = Taxonomy().read(self.taxonomy_file)
        
        fout = open(os.path.join(out_dir, prefix + '.classification.tsv'), 'w')   
        
        # We measure the mash distance for Genome placed on terminal branches
        analysed_nodes = []
        mash_dict = {}
        for nd in tree:
            list_subnode_initials = [subnd.taxon.__str__().replace("'",'')[0] for subnd in nd.leaf_iter()]
            list_subnode = [subnd.taxon.__str__().replace("'",'') for subnd in nd.leaf_iter()]
            if list_subnode_initials.count('U') == len(list_subnode_initials) - 1 and len(list_subnode_initials) > 1 and list_subnode[0] not in analysed_nodes:
                results = self._calculate_mash_distance(list_subnode,genomes)
                mash_dict = merge_two_dicts(mash_dict,results)
                analysed_nodes.extend(list_subnode)
        for k,v in mash_dict.iteritems():
            suffixed_name = add_ncbi_prefix(v.get("ref_genome"))
            taxa_str = ";".join(gtdb_taxonomy.get(suffixed_name))
            fout.write('%s\t%s\n' % (k, taxa_str))        
        
        user_genome_ids = set(read_fasta(user_msa_file).keys())
        user_genome_ids = user_genome_ids.difference(set(mash_dict.keys()))
        for leaf in tree.leaf_node_iter():
            if leaf.taxon.label in user_genome_ids:
                taxa = []
                cur_node = leaf
                while cur_node.parent_node:
                    _support, taxon, _aux_info = parse_label(cur_node.label)
                    if taxon:
                        for t in taxon.split(';')[::-1]:
                            taxa.append(t.strip())
                            
                    cur_node = cur_node.parent_node
                
                taxa_str = ';'.join(taxa[::-1])
                fout.write('%s\t%s\n' % (leaf.taxon.label, taxa_str))
        fout.close()
    
        
        print "THERE RESULTS SHOULD BE REFINED TO SEE IF A GENOME CAN BE ASSIGNED TO ANY SISTER TAXON"
        print "EX: PERHAPS THIS GENOME BELONGS TO A SISTER CLASS!"
        print "NEED TO GET FIXED PhyloRank THRESHOLDS INTO THE MIX."
        
        
    def _calculate_mash_distance(self,list_leaf,genomes):
        """ Calculate the Mash distance between all user genomes and the reference to classfy them at the species level"""
        try:
            self.tmp_output_dir = tempfile.mkdtemp()
            for leaf in list_leaf:
                if leaf.startswith('U'):
                    shutil.copy(genomes.get(leaf),self.tmp_output_dir)
            cmd = 'mash sketch -s 5000 -k 16 -o {0}/user_genomes {0}/*.fna -p {1} > /dev/null 2>&1'.format(self.tmp_output_dir,self.cpus)
            os.system(cmd)
            reference_db = os.path.join(Config.MASH_DIR,Config.MASH_DB)
            cmd = 'mash dist {0} {1}/user_genomes.msh -p {2} -d {3}> {1}/distances.tab'.format(reference_db,self.tmp_output_dir,self.cpus,Config.MASH_SPECIES_THRESHOLD) 
            os.system(cmd)
            if not os.path.isfile(os.path.join(self.tmp_output_dir,'user_genomes.msh')) or not os.path.isfile(os.path.join(self.tmp_output_dir,'distances.tab')):
                raise
            dict_parser_distance = self._parse_mash_results(os.path.join(self.tmp_output_dir,'distances.tab'))
            return dict_parser_distance
              
        except:
            if os.path.exists(self.tmp_output_dir):
                shutil.rmtree(self.tmp_output_dir)
            raise
        
    def _parse_mash_results(self,distance_file):
        dict_results = {}
        with open(distance_file) as distfile:
            for line in distfile:
                info = line.strip().split("\t")
                ref_genome = "_".join(info[0].split("_", 2)[:2])
                user_g = remove_extension(os.path.basename(info[1]))
                mash_dist = float(info[2])
                if user_g in dict_results:
                    if mash_dist < dict_results.get(user_g).get("mash_dist"):
                        dict_results[user_g]={"ref_genome":ref_genome,"mash_dist":mash_dist}
                else:
                    dict_results[user_g]={"ref_genome":ref_genome,"mash_dist":mash_dist}
                    
        return dict_results

                    
                
