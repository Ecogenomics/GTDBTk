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

from .HMMSequence import HMMSequence
from .HMMUnit import HMMUnit


class HMMResults(object):
    """
     This class has been adapted from the Perl module written by Genome Research Ltd.
     Perl authors: Rob Finn (rdf@sanger.ac.uk)
     Perl version: HMMResults.pm,v 1.3 2009-12-15 14:38:08 jt6 Exp $
     Python authors: Aaron Mussig (a.mussig@uq.edu.au)
     """

    def __init__(self):
        self.hmmerVersion = None
        self.hmmName = None
        self.seqDB = None
        self.hmmLength = None
        self.thisFile = None
        self.seedName = None
        self.seqs = dict()
        self.units = list()
        self.domThr = 25.0
        self.seqThr = 25.0
        self.evalueThr = None
        self.domTC = None
        self.seqTC = None
        self.domNC = None
        self.seqNC = None
        self.randSeedNum = None
        self.description = None
        self.seqName = None
        self.seqLength = None
        self.eof = False
        self.program = None

    def addHMMSeq(self, hmmSeq):
        """
        Adds a HMMSequence object to the results object
        """
        if not isinstance(hmmSeq, HMMSequence):
            sys.exit('Trying to add a non HMMSequence object')

        if hmmSeq.name in self.seqs:
            sys.exit('Trying to add the same sequence twice')

        self.seqs[hmmSeq.name] = hmmSeq

    def addHMMUnit(self, hmmUnit):
        """
        Adds HMM units (the actual region hit) to the HMMSequence in the object and for convenience to
        the results sets. All we store are duplicates of the references.
        """
        if not isinstance(hmmUnit, HMMUnit):
            sys.exit('Trying to add a non HMMUnit object')

        if self.seqs:
            if hmmUnit.name in self.seqs:
                self.seqs[hmmUnit.name].addHMMUnit(hmmUnit)
            else:
                sys.stderr.write('Could not add hmmUnit as the sequence has not been added\n')

        # More convenience we store the point to the hmmunit in an array
        self.units.append(hmmUnit)

    def remove_overlaps_by_clan(self, clanmap, nested):
        new = HMMResults()
        new.seqName = self.seqName

        for unit in sorted(self.units, key=lambda x: float(x.evalue), reverse=False):

            # check if it overlaps before adding
            o = None

            for u in new.units:
                if unit.name in clanmap and u.name in clanmap and clanmap[unit.name] == clanmap[u.name]:
                    if overlap(unit, u):
                        # TODO: Did i interpret the logic correctly here?
                        if unit.name in nested and u.name in nested[unit.name]:
                            continue
                        else:
                            o = True
                            break

            if not o:
                if unit.name not in new.seqs:
                    new_seq = HMMSequence()
                    new_seq.name = self.seqs[unit.name].name
                    new_seq.desc = self.seqs[unit.name].desc
                    new_seq.bits = self.seqs[unit.name].bits
                    new_seq.evalue = self.seqs[unit.name].evalue
                    new_seq.numberHits = self.seqs[unit.name].numberHits

                    new.addHMMSeq(new_seq)

                new.addHMMUnit(unit)

        return new


def overlap(unit1, unit2):
    """
    does unit1 overlap with unit2?
    """
    out = sorted([unit1, unit2], key=lambda x: int(x.seqFrom), reverse=False)
    u1 = out[0]
    u2 = out[1]

    if u2.seqFrom <= u1.seqTo:
        return True

    return False
