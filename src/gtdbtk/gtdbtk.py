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
import logging
import sys

from markers import Markers
from classify import Classify
from reroot_tree import RerootTree
import config.config as Config

from biolib.common import (check_dir_exists,
                            check_file_exists, 
                            make_sure_path_exists)
from biolib.taxonomy import Taxonomy
from biolib.external.execute import check_dependencies
from biolib.external.fasttree import FastTree 

class OptionsParser():

    def __init__(self, version):
        """Initialization."""
        
        self.version = version
        
        self.logger = logging.getLogger('timestamp')
        
    def _marker_set_id(self, bac120_ms, ar122_ms, rps23_ms):
        """Get unique identifier for marker set."""
        
        if bac120_ms:
            marker_set_id = "bac120"
        elif ar122_ms:
            marker_set_id = "ar122"
        elif rps23_ms:
            marker_set_id = "rps23"
            
        return marker_set_id
        
    def identify(self, options):
        """Identify marker genes in genomes."""
        
        try:
        
            if options.genome_dir:
                check_dir_exists(options.genome_dir)
                
            if options.batchfile:
                check_file_exists(options.batchfile)
                
            make_sure_path_exists(options.out_dir)
                
            markers = Markers(options.cpus)
            markers.identify(options.genome_dir,
                                options.batchfile,
                                options.proteins,
                                options.out_dir, 
                                options.prefix)
                                
            self.logger.info('Done.')
        
        except Exception as e:
            self.logger.info('GTDB-Tk has encountered an error.')
        
    def align(self, options):
        """Create MSA from marker genes."""
        
        if options.genome_dir:
            check_dir_exists(options.genome_dir)
            
        if options.batchfile:
            check_file_exists(options.batchfile)
            
        check_dir_exists(options.identify_dir)
        make_sure_path_exists(options.out_dir)
 
        marker_set_id = self._marker_set_id(options.bac120_ms,
                                            options.ar122_ms,
                                            options.rps23_ms)
          
        markers = Markers(options.threads)
        markers.align(options.genome_dir,
                        options.batchfile,
                        options.identify_dir,
                        marker_set_id,
                        options.taxa_filter,
                        options.min_perc_aa,
                        options.custom_msa_filters,
                        options.consensus,
                        options.min_perc_taxa,
                        options.out_dir,
                        options.prefix)
                       
        self.logger.info('Done.')
        
    def infer(self, options):
        """Infer tree from MSA."""
        
        check_file_exists(options.msa_file)
        make_sure_path_exists(options.out_dir)
        
        if (options.cpus > 1):
            check_dependencies(['FastTreeMP'])
        else:
            check_dependencies(['FastTree'])
        
        self.logger.info('Inferring tree with FastTree using %s+GAMMA.' % options.prot_model)
        fasttree = FastTree(multithreaded=(options.cpus > 1))

        tree_unrooted_output = os.path.join(options.out_dir, options.prefix + '.unrooted.tree')
        tree_log = os.path.join(options.out_dir, options.prefix + '.tree.log')
        tree_output_log = os.path.join(options.out_dir, 'fasttree.log')
        fasttree.run(options.msa_file, 
                        'prot', 
                        options.prot_model, 
                        tree_unrooted_output, 
                        tree_log, 
                        tree_output_log)
        
        self.logger.info('Done.')
        
    def classify(self, options):
        """Determine taxonomic classification of genomes."""
        
        if options.genome_dir:
            check_dir_exists(options.genome_dir)
            
        if options.batchfile:
            check_file_exists(options.batchfile)

        check_file_exists(options.user_msa_file)
        make_sure_path_exists(options.out_dir)
        
        marker_set_id = self._marker_set_id(options.bac120_ms,
                                            options.ar122_ms,
                                            options.rps23_ms)

        classify = Classify(options.cpus)
        classify.run(options.user_msa_file,
                     options.genome_dir,
                     options.batchfile,                    
                     marker_set_id,
                     options.out_dir,
                     options.prefix)
        
        self.logger.info('Done.')
        
    def root(self, options):
        """Root tree using outgroup."""
        
        check_file_exists(options.input_tree)
        
        gtdb_taxonomy = Taxonomy().read(Config.TAXONOMY_FILE)
        
        self.logger.info('Identifying genomes from the specified outgroup.')
        outgroup = set()
        for genome_id, taxa in gtdb_taxonomy.iteritems():
            if options.outgroup_taxon in taxa:
                outgroup.add(genome_id)

        reroot = RerootTree()
        reroot.root_with_outgroup(options.input_tree, 
                                    options.output_tree, 
                                    outgroup)
        
        self.logger.info('Done.')
        
    def decorate(self, options):
        """Decorate tree with GTDB taxonomy."""
        
        check_file_exists(options.input_tree)
        
        #Config.TAXONOMY_FILE
        self.logger.warning('NOT YET IMPLEMENTED!')
        
        self.logger.info('Done.')
        
    def parse_options(self, options):
        """Parse user options and call the correct pipeline(s)"""
        
        if(options.subparser_name == 'de_novo_wf'):
            check_dependencies(['prodigal', 'hmmalign'])
            if (options.cpus > 1):
                check_dependencies(['FastTreeMP'])
            else:
                check_dependencies(['FastTree'])
                
            self.identify(options)
            
            options.identify_dir = options.out_dir
            self.align(options)
            
            options.msa_file = os.path.join(options.out_dir, 
                                            options.prefix + ".msa.faa")
            self.infer(options)
            
            options.input_tree = os.path.join(options.out_dir, 
                                                options.prefix + ".unrooted.tree")
            options.output_tree = os.path.join(options.out_dir, 
                                                options.prefix + ".rooted.tree")
            self.root(options)
            
            self.decorate(options)
        elif(options.subparser_name == 'classify_wf'):
            check_dependencies(['prodigal', 'hmmalign', 'pplacer', 'guppy'])               
            self.identify(options)
            
            options.identify_dir = options.out_dir
            options.taxa_filter = None
            options.custom_msa_filters = False
            options.consensus = None
            options.min_perc_taxa = None
            self.align(options)
            
            options.user_msa_file = os.path.join(options.out_dir, 
                                                    options.prefix + ".user_msa.faa")
            self.classify(options)
        elif (options.subparser_name == 'identify'):
            self.identify(options)
        elif(options.subparser_name == 'align'):
            self.align(options)
        elif(options.subparser_name == 'infer'):
            self.infer(options)
        elif(options.subparser_name == 'classify'):
            self.classify(options)
        elif(options.subparser_name == 'root'):
            self.root(options)
        elif(options.subparser_name == 'decorate'):
            self.decorate(options)
        else:
            self.logger.error('Unknown GTDB-Tk command: "' + options.subparser_name + '"\n')
            sys.exit()
            
        return 0
