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

import sys

from .HMMMatch import HMMMatch
from .HMMUnit import HMMUnit


class HMMSequence(HMMMatch):
    """
    This class has been adapted from the Perl module written by Genome Research Ltd.
    Perl authors: Rob Finn (rdf@sanger.ac.uk), John Tate (jt6@sanger.ac.uk)
    Perl version: ?
    Python authors: Aaron Mussig (a.mussig@uq.edu.au)
    """

    def __init__(self):
        HMMMatch.__init__(self)
        self.sumEvalue = None
        self.H2mode = None
        self.sumScore = None
        self.desc = None
        self.numberHits = None
        self.exp = None
        self.hmmUnits = list()  # An array of HMMUnit

    def addHMMUnit(self, hmmUnit):
        """
        Adds a hmmUnit to a sequence. It checks that the variable passed in is a HMMUnit object
        """
        if isinstance(hmmUnit, HMMUnit):
            self.hmmUnits.append(hmmUnit)
        else:
            sys.stderr.write('%s is not a HMMUnit, not added\n.' % hmmUnit)
