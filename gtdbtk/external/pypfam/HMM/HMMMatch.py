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


class HMMMatch(object):
    """
    This class has been adapted from the Perl module written by Genome Research Ltd.
    Perl authors: Jaina Mistry (jm14@sanger.ac.uk), John Tate (jt6@sanger.ac.uk), Rob Finn (finnr@janelia.hhmi.org)
    Perl version: ?
    Python authors: Aaron Mussig (a.mussig@uq.edu.au)
    """

    def __init__(self):
        """
        evalue (float), bits (float), name (str) , bias (float)
        """
        self.evalue = None
        self.bits = None
        self.name = None
        self.bias = None
