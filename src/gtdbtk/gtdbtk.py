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

from markers import Markers
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
        
    def identify(self, options):
        """Identify marker genes in genomes."""
        
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
        
    def align(self, options):
        """Create MSA from marker genes."""
        
        if options.genome_dir:
            check_dir_exists(options.genome_dir)
            
        if options.batchfile:
            check_file_exists(options.batchfile)
            
        check_dir_exists(options.identify_dir)
        make_sure_path_exists(options.out_dir)
 
        marker_set_id = None
        if options.bac_ms:
            marker_set_id = "bac120"
        elif options.ar_ms:
            marker_set_id = "ar122"
          
        markers = Markers(options.threads)
        markers.align(options.genome_dir,
                        options.batchfile,
                        options.identify_dir,
                        marker_set_id,
                        options.taxa_filter,
                        options.min_perc_aa,
                        options.consensus,
                        options.min_perc_taxa,
                        options.out_dir,
                        options.prefix)
                       
        self.logger.info('Done.')
        
    def infer(self, options):
        """Infer tree from MSA."""
        
        check_file_exists(options.msa_file)
        
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
        
        if (options.subparser_name == 'identify'):
            self.identify(options)
        elif(options.subparser_name == 'align'):
            self.align(options)
        elif(options.subparser_name == 'infer'):
            self.infer(options)
        elif(options.subparser_name == 'root'):
            self.root(options)
        elif(options.subparser_name == 'decorate'):
            self.decorate(options)
        else:
            self.logger.error('Unknown GTDB-Tk command: "' + args.subparser_name + '"\n')
            sys.exit()
            
        return 0
