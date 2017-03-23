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

from TreeManager import TreeManager
from markers import Markers

from biolib.common import (check_dir_exists,
                            check_file_exists, 
                            make_sure_path_exists)


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
        """Create multiple sequence alignment from marker genes"""
        
        check_file_exists(options.batchfile)
        make_sure_path_exists(options.in_dir)
 
        domain = None
        if options.bac_domain:
            domain = "bacteria"
        elif options.arc_domain:
            domain = "archaea"
          
        markers = Markers(options.threads)
        markers.align(options.batchfile,
                       options.in_dir,
                       domain,
                       options.filter_taxa,
                       options.min_perc_aa,
                       options.consensus,
                       options.min_perc_taxa,
                       options.out_dir,
                       options.prefix)
                       
        self.logger.info('Done.')
        
    def parse_options(self, options):
        """Parse user options and call the correct pipeline(s)"""
        
        if (options.subparser_name == 'identify'):
            self.identify(options)
        elif(options.subparser_name == 'align'):
            self.align(options)
        else:
            self.logger.error('Unknown GTDB-Tk command: "' + args.subparser_name + '"\n')
            sys.exit()
            
        return 0
