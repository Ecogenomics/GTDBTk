#!/usr/bin/env python
###############################################################################
#                                                                             #
#                                                                             #
#                                                                             #
#    Entry point. See gtdbtk/main.py for internals                             #
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

from TreeManager import TreeManager
from GtdbManager import GtdbManager


class GtdbTKOptionsParser():

    def __init__(self, version):
        self.GTVersion = version

    def parseOptions(self, options):
        if(options.subparser_name == 'align'):
            # parse raw input
            print "*******************************************************************************"
            print " [[GtdbTK %s]] Aligning Genomes..." % self.GTVersion
            print "*******************************************************************************"
            gtdb_mngr = GtdbManager(options.threads)
            domain = None
            if options.bac_domain:
                domain = "bacteria"
            if options.arc_domain:
                domain = "archaea"
            success = gtdb_mngr.AlignedGenomes(options.batchfile,
                                               options.indir,
                                               domain,
                                               options.filter_taxa,
                                               options.min_perc_aa,
                                               options.consensus,
                                               options.min_perc_taxa,
                                               options.out_dir,
                                               options.prefix)
            if not success:
                print "ERROR"

        elif (options.subparser_name == 'identify'):
            # parse raw input
            print "*******************************************************************************"
            print " [[GtdbTK %s]] Identifying genome domains..." % self.GTVersion
            print "*******************************************************************************"
            gtdb_mngr = GtdbManager(options.threads)
            success = gtdb_mngr.IdentifyMarkers(options.batchfile, options.outdir, options.prefix)
            if not success:
                print "ERROR"

        return 0
