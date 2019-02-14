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

from .HMMMatch import HMMMatch


class HMMUnit(HMMMatch):
    """
    This class has been adapted from the Perl module written by Genome Research Ltd.
    Perl authors: Rob Finn (rdf@sanger.ac.uk), John Tate (jt6@sanger.ac.uk)
    Perl version: ?
    Python authors: Aaron Mussig (a.mussig@uq.edu.au)
    """

    def __init__(self):
        HMMMatch.__init__(self)
        self.proteinCoos = None
        self.seqEvalue = None
        self.domain = None
        self.seqFrom = None
        self.seqTo = None
        self.domEvalue = None
        self.hmmalign = dict()
        self.hmmFrom = None
        self.hmmTo = None
        self.envFrom = None
        self.envTo = None
        self.coreFrom = None
        self.coreTo = None
        self.aliAcc = None
        self.sig = None
