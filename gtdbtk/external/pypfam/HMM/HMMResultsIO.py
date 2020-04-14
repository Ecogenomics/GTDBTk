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

import re
import sys
from collections import deque

from .HMMResults import HMMResults
from .HMMSequence import HMMSequence
from .HMMUnit import HMMUnit


class HMMResultsIO(object):
    """
    This class has been adapted from the Perl module written by Genome Research Ltd.
    Perl authors: Rob Finn (rdf@sanger.ac.uk), John Tate (jt6@sanger.ac.uk)
    Perl version: $Id: HMMResultsIO.pm,v 1.2 2009-12-01 15:42:20 jt6 Exp $
    Python authors: Aaron Mussig (a.mussig@uq.edu.au)
    """

    def __init__(self):
        self.align = 0
        self.outfile = 'OUTPUT'
        self.pfamout = 'PFAMOUT'
        self.scores = 'scores'

        # Python adaptation: Pre-compile regex patterns used in _readHeader
        self.re_header_end = re.compile(r'^Scores for complete')
        self.re_hmm_name = re.compile(r'^# query HMM file:\s+(\S+)')
        self.re_seq_db = re.compile(r'^# target sequence database:\s+(\S+)')
        self.re_this_file = re.compile(r'^output directed to file:\s+(\S+)')
        self.re_seed_name__hmm_length = re.compile(r'^Query:\s+(\S+)\s+\[M\=(\d+)\]')
        self.re_seq_name__seq_length = re.compile(r'^Query:\s+(\S+)\s+\[L\=(\d+)\]')
        self.re_e_value_Thr = re.compile(r'^sequence E-value threshold: <= (\d+)')
        self.re_rand_seed_num = re.compile(r'^# Random generator seed:      (\d+)')
        self.re_descr = re.compile(r'^Description:\s+(.*)')
        self.re_program = re.compile(r'^# (phmmer|hmmsearch|hmmscan|jackhmmer)')
        self.re_header_skip = re.compile(r'(^#)|(^$)')
        self.re_header_skip_2 = re.compile(r'^Accession')
        self.re_eof = re.compile(r'^\[ok\]')

        # Python adaptation: Pre-compile the regex used in parsing the sequence hits.
        self.re_seq_hit_1 = re.compile(r'^Domain annotation for each [sequence|model]')
        self.re_seq_hit_2 = re.compile(r'^Domain and alignment annotation for each [sequence|model]')
        self.re_seq_hit_3 = re.compile(r'^\s+(E-value|---)')
        self.re_seq_hit_4 = re.compile(r'^$')
        self.re_seq_hit_5 = re.compile(r'No hits detected that satisfy reporting thresholds')
        self.re_seq_hit_6 = re.compile(r'\s+')

        # Python adaptation: Pre-compile the regex used in parsing the domain hits.
        self.re_dom_1 = re.compile(r'^Internal')
        self.re_dom_2 = re.compile(r'\>\>\s+(\S+)')

        # Python adaptation: Pre-compile the regex used in parsing the unit data.
        self.re_unitd_1 = re.compile(r'^[(\/\/|Internal)]')
        self.re_unitd_2 = re.compile(r'^\>\>\s+(\S+)')
        self.re_unitd_3 = re.compile(r'^\s+Alignments for each domain:')
        self.re_unitd_4 = re.compile(r'^\s+(#\s+score|---)')
        self.re_unitd_5 = re.compile(r'^$')
        self.re_unitd_6 = re.compile(r'^\s+\d+\s+')
        self.re_unitd_7 = re.compile(r'\s+')
        self.re_unitd_8 = re.compile(r'^\s+\[No individual domains')
        self.re_unitd_9 = re.compile(r'^\s+([x\.]+)\s+RF$')
        self.re_unitd_10 = re.compile(r'^\s+([0-9\*\.]+)\s+PP$')
        self.re_unitd_11 = re.compile(r'^\s+(\S+)\s+CS$')
        self.re_unitd_12 = re.compile(r'^\s+==\s+domain\s+(\d+)')
        self.re_unitd_13 = re.compile(r'^\s+(.*)\s*$')  # Modified from Perl version as no newline break in iterator.
        self.re_unitd_14 = re.compile(r'^$')
        self.re_unitd_15 = re.compile(r'^[(\/\/|Internal)]')
        self.re_unitd_16 = re.compile(r'^\>\>\s+(\S+)')

        # Python adaptation: Pre-compile the regex used in parsing the footer.
        self.re_footer_1 = re.compile(r'\/\/')

        # Python adaptation: Pre-compile the regex used in the ascii_output function
        self._re_ascii_1 = re.compile(r'Pfam\-B')

    def parseMultiHMMER3(self, filename):

        fh = deque(filename.splitlines())

        hmmResAll = list()
        program = None

        while len(fh) > 0:
            hmmRes = HMMResults()

            eof = self._readHeader(fh, hmmRes)
            if eof:
                break

            hmmResAll.append(hmmRes)

            if hmmRes.program:
                program = hmmRes.program
            else:
                hmmRes.program = program

            self._readSeqHits(fh, hmmRes)
            self._readUnitHits(fh, hmmRes)
            self._readFooter(fh, hmmRes)

        return hmmResAll

    def _readHeader(self, fh, hmmRes):
        """
        Reads the header section from a HMMER3 hmmsearch
        """

        while len(fh) > 0:
            hs = fh.popleft()

            # Python adaptation: Execute all regex searches which extract information.
            hmm_name = self.re_hmm_name.search(hs)
            seq_db = self.re_seq_db.search(hs)
            this_file = self.re_this_file.search(hs)
            seed_name__hmm_length = self.re_seed_name__hmm_length.search(hs)
            seq_name__seq_length = self.re_seq_name__seq_length.search(hs)
            e_value_Thr = self.re_e_value_Thr.search(hs)
            rand_seed_num = self.re_rand_seed_num.search(hs)
            descr = self.re_descr.search(hs)
            program = self.re_program.search(hs)

            if self.re_header_end.search(hs):
                break
            elif hmm_name:
                hmmRes.hmmName = hmm_name.group(1)
            elif seq_db:
                hmmRes.seqDB = seq_db.group(1)
            elif this_file:
                hmmRes.thisFile = this_file.group(1)
            elif seed_name__hmm_length:
                hmmRes.seedName = seed_name__hmm_length.group(1)
                hmmRes.hmmLength = seed_name__hmm_length.group(2)
            elif seq_name__seq_length:
                hmmRes.seqName = seq_name__seq_length.group(1)
                hmmRes.seqLength = seq_name__seq_length.group(2)
            elif e_value_Thr:
                hmmRes.evalueThr = e_value_Thr.group(1)
            elif rand_seed_num:
                hmmRes.randSeedNum = rand_seed_num.group(1)
            elif descr:
                hmmRes.description = descr.group(1)
            elif program:
                hmmRes.program = program.group(1)
            elif self.re_header_skip.search(hs):
                continue
            elif self.re_header_skip_2.search(hs):
                continue
            elif self.re_eof.search(hs):
                return True
            else:
                sys.exit("Failed to parse hmmsearch results in header section\n")

    def _readSeqHits(self, fh, hmmRes):
        """
        Reads the sequence hits from a HMMER3 hmmsearch
        """

        while len(fh) > 0:
            hs = fh.popleft()

            if self.re_seq_hit_1.search(hs):
                break
            elif self.re_seq_hit_2.search(hs):
                break
            elif self.re_seq_hit_3.search(hs):
                continue
            elif self.re_seq_hit_4.search(hs):
                continue
            elif self.re_seq_hit_5.search(hs):
                continue
            else:

                # Assume that we have a sequence match
                re_6 = self.re_seq_hit_6.split(hs)

                # CHeck if it is what is expected
                if len(re_6) < 10:
                    sys.exit('Expected at least 10 pieces of data.\n')

                descr = ' '.join(re_6[10:])

                hmmSeq = HMMSequence()
                hmmSeq.evalue = float(re_6[1])
                hmmSeq.bits = float(re_6[2])
                hmmSeq.bias = float(re_6[3])
                hmmSeq.exp = float(re_6[7])
                hmmSeq.numberHits = int(re_6[8])
                hmmSeq.name = re_6[9]
                hmmSeq.desc = descr if descr else '-'
                hmmRes.addHMMSeq(hmmSeq)

    def _readUnitHits(self, fh, hmmRes):
        """
        Reads the unit (domain) hits from a HMMER3 hmmsearch
        """

        if hmmRes.eof:
            return

        # Parse the domain hits section
        while len(fh) > 0:
            hs = fh.popleft()

            # Perform the regex searches
            re_dom_2 = self.re_dom_2.search(hs)

            if self.re_dom_1.search(hs):
                break

            elif re_dom_2:
                seqId = re_dom_2.group(1)
                self._readUnitData(seqId, fh, hmmRes)

                if hmmRes.eof:
                    return

    def _readUnitData(self, seqId, fh, hmmRes):

        if hmmRes.eof:
            return

        hmmName = hmmRes.seedName
        seqName = hmmRes.seqName

        units = list()
        align = True
        recurse = False
        eof = False
        nextSeqId = None

        # Parse the domain hits section
        while len(fh) > 0:
            hs = fh.popleft()

            # Run the regex searches which generate output
            unitd_2 = self.re_unitd_2.search(hs)

            if self.re_unitd_1.search(hs):
                align = False
                recurse = False
                eof = True
                break

            elif unitd_2:
                nextSeqId = unitd_2.group(1)
                align = False
                recurse = True
                break

            elif self.re_unitd_3.search(hs):
                align = True
                recurse = False
                break

            elif self.re_unitd_4.search(hs):
                # Two human readable lines
                continue

            elif self.re_unitd_5.search(hs):
                # blank line
                continue

            elif self.re_unitd_6.search(hs):
                dMatch = self.re_unitd_7.split(hs)
                if len(dMatch) != 17:
                    sys.exit('Expected 16 elements of data.')

                hmmUnit = HMMUnit()
                hmmUnit.name = seqId
                hmmUnit.domain = dMatch[1]
                hmmUnit.bits = float(dMatch[3])
                hmmUnit.bias = float(dMatch[4])
                hmmUnit.domEvalue = float(dMatch[5])
                hmmUnit.evalue = float(dMatch[6])
                hmmUnit.hmmFrom = int(dMatch[7])
                hmmUnit.hmmTo = int(dMatch[8])
                hmmUnit.seqFrom = int(dMatch[10])
                hmmUnit.seqTo = int(dMatch[11])
                hmmUnit.envFrom = int(dMatch[13])
                hmmUnit.envTo = int(dMatch[14])
                hmmUnit.aliAcc = float(dMatch[16])

                units.append(hmmUnit)
                continue

            elif self.re_unitd_8.search(hs):
                align = False
                continue

            else:
                sys.exit('Did not parse line %s' % hs)

        '''
        #  == domain 1    score: 244.0 bits;  conditional E-value: 9.5e-76
        #      SEED   1 medrlallkaisasakdlvalaasrGaksipspvkttavkfdplptPdldalrtrlkeaklPakaiksalsayekaCarWrsdleeafdktaksvsPanlhllealrirlyteqvekWlvqvlevaerWkaemekqrahiaatmgp 146
        #               m+++la+l++isa+akd++ala+srGa+++ +p++tt+++fd+l++P+ld++rtrl+ea+lP+kaik++lsaye+aCarW++dleeafd+ta+s+sP+n+++l++lr+rly+eqv+kWl++vl+v+erWkaemekqrahi+atmgp
        #  P37935.1   1 MAELLACLQSISAHAKDMMALARSRGATGS-RPTPTTLPHFDELLPPNLDFVRTRLQEARLPPKAIKGTLSAYESACARWKHDLEEAFDRTAHSISPHNFQRLAQLRTRLYVEQVQKWLYEVLQVPERWKAEMEKQRAHINATMGP 145
        #               899***************************.******************************************************************************************************************8 PP
        #
        # OR....
        #
        #  == domain 1    score: 27.6 bits;  conditional E-value: 7.4e-10
        #   PF00018  17 LsfkkGdvitvleksee.eWwkaelkdg.keGlvPsnYvep 55 
        #               L++++Gd+++++++++e++Ww++++++++++G++P+n+v+p
        #  P15498.4 617 LRLNPGDIVELTKAEAEqNWWEGRNTSTnEIGWFPCNRVKP 657
        #               7899**********9999*******************9987 PP
        '''

        if align:

            # Specifically for python
            pattern1 = None
            pattern2 = None

            if hmmName and hmmRes.program == 'hmmsearch':

                pattern1 = re.compile(r'^\s+%s\s+\d+\s+(\S+)\s+\d+' % hmmName)
                seqId = re.sub('(\W)', r'\\\1', seqId)
                # $id =~ s/\|/\\|/g;  #Escape '|', '[' and ']' characters
                # $id =~ s/\[/\\[/g;
                # $id =~ s/\]/\\]/g;
                pattern2 = re.compile(r'^\s+%s\s+\d+\s+(\S+)\s+\d+' % seqId)

            elif seqName and hmmRes.program == 'hmmscan':
                tmpSeqName = seqName
                tmpSeqName = re.sub('(\W)', r'\\\1', tmpSeqName)

                pattern1 = re.compile(r'^\s+%s\s+\d+\s+(\S+)\s+\d+' % seqId)
                pattern2 = re.compile(r'^\s+%s\s+\d+\s+(\S+)\s+\d+' % tmpSeqName)

            elif seqName and (hmmRes.program == 'phmmer' or hmmRes.program == 'jackhmmer'):
                sys.exit("seqName and (hmmRes.program == 'phmmer' or hmmRes.program == 'jackhmmer' is not implemented.")

            recurse = False
            matchNo = None
            hmmlen = 0

            while len(fh) > 0:
                hs = fh.popleft()

                # Run a search for each of the patterns.
                pattern1_res = pattern1.search(hs)
                pattern2_res = pattern2.search(hs)
                re_unitd_9_res = self.re_unitd_9.search(hs)
                re_unitd_10_res = self.re_unitd_10.search(hs)
                re_unitd_11_res = self.re_unitd_11.search(hs)
                re_unitd_12_res = self.re_unitd_12.search(hs)
                re_unitd_13_res = self.re_unitd_13.search(hs)
                re_unitd_14_res = self.re_unitd_14.search(hs)
                re_unitd_15_res = self.re_unitd_15.search(hs)
                re_unitd_16_res = self.re_unitd_16.search(hs)

                if pattern1_res:
                    dict_hmmalign = units[matchNo - 1].hmmalign
                    if 'hmm' in dict_hmmalign:
                        dict_hmmalign['hmm'] += pattern1_res.group(1)
                    else:
                        dict_hmmalign['hmm'] = pattern1_res.group(1)

                    hmmlen = len(pattern1_res.group(1))

                elif pattern2_res:
                    dict_hmmalign = units[matchNo - 1].hmmalign
                    if 'seq' in dict_hmmalign:
                        dict_hmmalign['seq'] += pattern2_res.group(1)
                    else:
                        dict_hmmalign['seq'] = pattern2_res.group(1)

                # ^\s+([x\.]+)\s+RF$
                elif re_unitd_9_res:
                    rf = re_unitd_9_res.group(1)
                    dict_hmmalign = units[matchNo - 1].hmmalign
                    if 'rf' in dict_hmmalign:
                        dict_hmmalign['rf'] += rf
                    else:
                        dict_hmmalign['rf'] = rf

                # ^\s+([0-9\*\.]+)\s+PP$
                elif re_unitd_10_res:
                    pp = re_unitd_10_res.group(1)

                    dict_hmmalign = units[matchNo - 1].hmmalign
                    if 'pp' in dict_hmmalign:
                        dict_hmmalign['pp'] += pp
                    else:
                        dict_hmmalign['pp'] = pp

                # ^\s+(\S+)\s+CS$
                elif re_unitd_11_res:
                    cs = re_unitd_11_res.group(1)
                    dict_hmmalign = units[matchNo - 1].hmmalign
                    if 'cs' in dict_hmmalign:
                        dict_hmmalign['cs'] += cs
                    else:
                        dict_hmmalign['cs'] = cs

                # ^\s+==\s+domain\s+(\d+)
                elif re_unitd_12_res:
                    matchNo = int(re_unitd_12_res.group(1))

                # ^\s+(.*)\s+$
                elif re_unitd_13_res:
                    hs = hs.rstrip()
                    m1 = hs[-hmmlen:]

                # ^$
                elif re_unitd_14_res:
                    continue

                # ^[(\/\/|Internal)]
                elif re_unitd_15_res:
                    align = False
                    recurse = False
                    eof = True
                    break

                # ^\>\>\s+(\S+)
                elif re_unitd_16_res:
                    nextSeqId = re_unitd_16_res.group(1)
                    recurse = True
                    break

                else:
                    sys.exit('Did not parse %s in units' % hs)

        # foreach my u (@units)
        for u in units:
            hmmRes.addHMMUnit(u)

        hmmRes.eof = eof

        if recurse and nextSeqId:
            self._readUnitData(nextSeqId, fh, hmmRes)

        return

    def _readFooter(self, fh, hmmRes):

        # We are going to parse something like this!

        #  Internal pipeline statistics summary:
        # -------------------------------------
        # Query sequence(s):                         1  (360 residues)
        # Target model(s):                           7  (836 nodes)
        # Passed MSV filter:                         2  (0.285714); expected 0.1 (0.02)
        # Passed Vit filter:                         1  (0.142857); expected 0.0 (0.001)
        # Passed Fwd filter:                         1  (0.142857); expected 0.0 (1e-05)
        # Initial search space (Z):                  7  [actual number of targets]
        # Domain search space  (domZ):               1  [number of targets reported over threshold]
        ## CPU time: 0.00u 0.00s 00:00:00.00 Elapsed: 00:00:00
        ## Mc/sec: inf
        # //

        while len(fh) > 0:
            hs = fh.popleft()

            if self.re_footer_1.search(hs):
                break

    def write_ascii_out(self, HMMResults_in, fh, scanData, e_seq, e_dom, b_seq, b_dom):

        if not scanData._max_seqname or scanData._max_seqname < 1:
            scanData._max_seqname = 20

        ga = None
        if e_seq or e_dom:
            sys.exit('FATAL: e_seq/e_dom functionality has not been implemented in Python')
        elif b_seq or b_dom:
            sys.exit('FATAL: b_seq/b_dom functionality has not been implemented in Python')
        else:
            ga = True

        for unit in sorted(HMMResults_in.units, key=lambda x: x.seqFrom, reverse=False):

            # Pfam\-B
            if self.re_footer_1.search(unit.name):
                sys.exit('FATAL: Pfam-B functionality has not been implemented in Python')

            else:

                # Filter results based on thresholds
                if ga:
                    if not unit.sig:
                        continue
                if e_seq:
                    sys.exit('FATAL: e_seq functionality has not been implemented in Python')
                if b_seq:
                    sys.exit('FATAL: b_seq functionality has not been implemented in Python')

                clan = "No_clan" if unit.name not in scanData._clanmap else scanData._clanmap[unit.name]

                seq_name_padded = HMMResults_in.seqName + ' ' * (scanData._max_seqname - len(HMMResults_in.seqName))
                fh.write('%s %6d %6d %6d %6d %-11s %-16s %7s %5d %5d %5d %8s %9s %3d %-8s ' %
                         (seq_name_padded,
                          unit.seqFrom,
                          unit.seqTo,
                          unit.envFrom,
                          unit.envTo,
                          scanData._accmap[unit.name],
                          unit.name,
                          scanData._type[unit.name],
                          unit.hmmFrom,
                          unit.hmmTo,
                          scanData._model_len[unit.name],
                          unit.bits,
                          unit.evalue,
                          unit.sig,
                          clan))

                fh.write('\n')

                if scanData._translate:
                    sys.exit('FATAL: translate functionality has not been implemented in Python')

            if scanData._align:
                sys.exit('FATAL: Align functionality has not been implemented in Python')
