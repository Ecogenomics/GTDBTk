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
import re
import subprocess
import sys
from datetime import datetime

from ..HMM.HMMResults import HMMResults
from ..HMM.HMMResultsIO import HMMResultsIO
from ..HMM.HMMSequence import HMMSequence
from ..HMM.HMMUnit import HMMUnit


class PfamScan(object):
    """
    This class has been adapted from the Perl module written by Genome Research Ltd.
    Perl authors: Jaina Mistry (jm14@sanger.ac.uk), John Tate (jt6@sanger.ac.uk), Rob Finn (finnr@janelia.hhmi.org)
    Perl version: PfamScan.pm,v 1.4 2010-01-12 09:41:42 jm14 Exp $
    Python authors: Aaron Mussig (a.mussig@uq.edu.au)
    """

    def __init__(self, **kwargs):
        """
        The only constructor for the object. Accepts a set of arguments that specify the parameters for the search:
        cut_off, cir, clan_overlap, fasta, sequence, align, hmm, as
        """

        # Python adaptation: Initialise variables
        # Function: _read_pfam_data
        self._read = dict()

        # To avoid hard coding the location for the binary, we assume it will be on the path.....
        self._HMMSCAN = 'hmmscan'

        # handle arguments, if we were given any here
        self._process_args(kwargs)

    def _process_args(self, kwargs):
        """
        Handles the input arguments.
        """

        # make sure we get a sequence
        if 'fasta' in kwargs and 'sequence' in kwargs:
            sys.exit('FATAL: "-fasta" and "-sequence" are mutually exclusive')
        elif 'fasta' in kwargs:
            if not os.path.isfile(kwargs['fasta']):
                sys.exit('FATAL: fasta file "%s" doesn\'t exist' % kwargs['fasta'])
        elif 'sequence' in kwargs:
            if len(kwargs['sequence'] < 1):
                sys.exit('FATAL: no sequence given')
        else:
            sys.exit('FATAL: must specify either "-fasta" or "-sequence"')

        # check the cut off
        if ('e_seq' in kwargs and ('b_seq' in kwargs or 'b_dom' in kwargs)) or \
                ('b_seq' in kwargs and ('e_seq' in kwargs or 'e_dom' in kwargs)) or \
                ('b_dom' in kwargs and 'e_dom' in kwargs):
            sys.exit('FATAL: can\'t use e value and bit score threshold together')

        self._hmmscan_cutoff = list()
        if 'e_seq' in kwargs:
            sys.exit('FATAL: -e_seq functionality has not been implemented in Python')
        if 'e_dom' in kwargs:
            sys.exit('FATAL: -e_dom functionality has not been implemented in Python')
        if 'b_seq' in kwargs:
            sys.exit('FATAL: -b_seq functionality has not been implemented in Python')
        if 'b_dom' in kwargs:
            sys.exit('FATAL: -b_dom functionality has not been implemented in Python')

        if not self._hmmscan_cutoff:
            self._hmmscan_cutoff.append('--cut_ga')

        # make sure we have a valid directory for the HMM data files
        if not os.path.isdir(kwargs['dir']):
            sys.exit('FATAL: directory "%s" does not exist' % kwargs['dir'])

        # Populate the object
        self._cut_off = kwargs.get('cut_off', None)
        self._dir = kwargs.get('dir', None)
        self._clan_overlap = kwargs.get('clan_overlap', None)
        self._fasta = kwargs.get('fasta', None)
        self._align = kwargs.get('align', None)
        self._as = kwargs.get('as', None)
        self._sequence = kwargs.get('sequence', None)
        self._cpu = kwargs.get('cpu', None)
        self._translate = kwargs.get('translate', None)

        self._hmmlib = list()
        if 'hmmlib' in kwargs:
            if type(kwargs['hmmlib']) == list:
                [self._hmmlib.append(x) for x in kwargs['hmmlib']]
            else:
                self._hmmlib.append(kwargs['hmmlib'])
        else:
            self._hmmlib.append('Pfam-A.hmm')

        # Now check that the library exists in the data dir!
        for hmmlib in self._hmmlib:

            if not os.path.isfile(os.path.join(self._dir, hmmlib)) and \
                    os.path.isfile(os.path.join(self._dir, hmmlib + '.h3f')) and \
                    os.path.isfile(os.path.join(self._dir, hmmlib + '.h3i')) and \
                    os.path.isfile(os.path.join(self._dir, hmmlib + '.h3m')) and \
                    os.path.isfile(os.path.join(self._dir, hmmlib + '.h3p')) and \
                    os.path.isfile(os.path.join(self._dir, hmmlib + '.dat')):
                sys.exit('FATAL: can\'t find %s and/or %s binaries in "%s' % (hmmlib, hmmlib, kwargs['dir']))

            # read the necessary data, if it's not been read already
            self._read_pfam_data()

        self._max_seqname = 0

        # if there's nothing in "_sequence" try to load a fasta file
        if not self._sequence:
            self._read_fasta()

        # check again for a sequence. If we don't have one at this point, bail with
        # an error
        if not self._sequence:
            sys.exit('FATAL: no sequence given')

        # read fasta file, store maximum sequence name and store sequences for active
        # sites prediction
        if not self._max_seqname:
            self._parse_sequence()

        if self._as:
            sys.exit('FATAL: -as functionality has not been implemented in Python')

        if self._translate:
            sys.exit('FATAL: -translate functionality has not been implemented in Python')

        # see if a version number was specified
        self._version = kwargs.get('version', None)

    def _read_pfam_data(self):
        """
        Reads the Pfam data file ("Pfam-A.scan.dat") and populates the C<accmap>,
        C<nested> and C<clanmap> hashes on the object.
        """

        # print STDERR "reading " . $self->{_hmmlib} . ".dat\n" if($ENV{DEBUG});
        self._accmap = dict()
        self._nested = dict()
        self._clanmap = dict()
        self._desc = dict()
        self._seqGA = dict()
        self._domGA = dict()
        self._type = dict()
        self._model_len = dict()

        # Python adaptation: Pre-compile regex
        re_read_pfam_1 = re.compile(r'^\#=GF ID\s+(\S+)')
        re_read_pfam_2 = re.compile(r'^\#=GF\s+AC\s+(\S+)')
        re_read_pfam_3 = re.compile(r'^\#=GF\s+DE\s+(.+)')
        re_read_pfam_4 = re.compile(r'^\#=GF\s+GA\s+(\S+)\;\s+(\S+)\;')
        re_read_pfam_5 = re.compile(r'^\#=GF\s+TP\s+(\S+)')
        re_read_pfam_6 = re.compile(r'^\#=GF\s+ML\s+(\d+)')
        re_read_pfam_7 = re.compile(r'^\#=GF\s+NE\s+(\S+)')
        re_read_pfam_8 = re.compile(r'^\#=GF\s+CL\s+(\S+)')

        for hmmlib in self._hmmlib:
            scandat = '%s/%s.dat' % (self._dir, hmmlib)

            try:
                with open(scandat, 'r') as f:
                    SCANDAT = f.readlines()
            except IOError as e:
                sys.exit('FATAL: Couldn\'t open "%s" data file: %s' % (scandat, e.message))

            v_id = None
            for line in SCANDAT:

                # Python adaptation: Pre-search all regex patterns, as can't search and return in-place.
                res_read_pfam_1 = re_read_pfam_1.search(line)  # ^\#=GF ID\s+(\S+)
                res_read_pfam_2 = re_read_pfam_2.search(line)  # ^\#=GF\s+AC\s+(\S+)
                res_read_pfam_3 = re_read_pfam_3.search(line)  # ^\#=GF\s+DE\s+(.+)
                res_read_pfam_4 = re_read_pfam_4.search(line)  # ^\#=GF\s+GA\s+(\S+)\;\s+(\S+)\;
                res_read_pfam_5 = re_read_pfam_5.search(line)  # ^\#=GF\s+TP\s+(\S+)
                res_read_pfam_6 = re_read_pfam_6.search(line)  # ^\#=GF\s+ML\s+(\d+)
                res_read_pfam_7 = re_read_pfam_7.search(line)  # ^\#=GF\s+NE\s+(\S+)
                res_read_pfam_8 = re_read_pfam_8.search(line)  # ^\#=GF\s+CL\s+(\S+)

                if res_read_pfam_1:
                    v_id = res_read_pfam_1.group(1)
                elif res_read_pfam_2:
                    self._accmap[v_id] = res_read_pfam_2.group(1)
                elif res_read_pfam_3:
                    self._desc[v_id] = res_read_pfam_3.group(1)
                elif res_read_pfam_4:
                    self._seqGA[v_id] = float(res_read_pfam_4.group(1))
                    self._domGA[v_id] = float(res_read_pfam_4.group(2))
                elif res_read_pfam_5:
                    self._type[v_id] = res_read_pfam_5.group(1)
                elif res_read_pfam_6:
                    self._model_len[v_id] = int(res_read_pfam_6.group(1))
                elif res_read_pfam_7:
                    self._nested[v_id] = {res_read_pfam_7.group(1): 1}
                    self._nested[res_read_pfam_7.group(1)] = {v_id: 1}
                elif res_read_pfam_8:
                    self._clanmap[v_id] = res_read_pfam_8.group(1)

            # set a flag to show that we've read the data files already
            self._read[hmmlib] = True

    def _convert_results_search_to_scan(self, search_results):
        """
        Converts the search format to the scan format
        """

        scan_results = dict()

        for search_result in search_results:
            for seq_id in search_result.seqs:
                this_seq_obj = search_result.seqs[seq_id]

                if seq_id not in scan_results:
                    new_scan_result = HMMResults()
                    new_scan_result.seqName = this_seq_obj.name
                    new_scan_result.description = this_seq_obj.desc
                    new_scan_result.program = search_result.program

                    scan_results[seq_id] = new_scan_result

                this_scan_result = scan_results[seq_id]

                new_hmm_seq = HMMSequence()
                new_hmm_seq.evalue = this_seq_obj.evalue
                new_hmm_seq.bits = this_seq_obj.bits
                new_hmm_seq.bias = this_seq_obj.bias
                new_hmm_seq.exp = this_seq_obj.exp
                new_hmm_seq.numberHits = this_seq_obj.numberHits
                new_hmm_seq.name = search_result.seedName
                new_hmm_seq.desc = search_result.description
                this_scan_result.addHMMSeq(new_hmm_seq)

                for search_unit in this_seq_obj.hmmUnits:
                    new_hmm_unit = HMMUnit()
                    new_hmm_unit.name = search_result.seedName
                    new_hmm_unit.domain = search_unit.domain
                    new_hmm_unit.hmmalign = search_unit.hmmalign
                    new_hmm_unit.bits = search_unit.bits
                    new_hmm_unit.bias = search_unit.bias
                    new_hmm_unit.domEvalue = search_unit.domEvalue
                    new_hmm_unit.evalue = search_unit.evalue
                    new_hmm_unit.hmmFrom = search_unit.hmmFrom
                    new_hmm_unit.hmmTo = search_unit.hmmTo
                    new_hmm_unit.seqFrom = search_unit.seqFrom
                    new_hmm_unit.seqTo = search_unit.seqTo
                    new_hmm_unit.envFrom = search_unit.envFrom
                    new_hmm_unit.envTo = search_unit.envTo
                    new_hmm_unit.aliAcc = search_unit.aliAcc
                    this_scan_result.addHMMUnit(new_hmm_unit)

                this_scan_result.eof = True

        values = list()
        for key, value in sorted(scan_results.items()):
            values.append(value)

        return values

    def _read_fasta(self):
        """
        Reads a sequence from the fasta-format file that was specified in the parameters.
        """
        try:
            with open(self._fasta, 'r') as f:
                rows = f.readlines()
        except IOError as e:
            sys.exit('FATAL: Couldn\'t open fasta file "%s" %s' % (self._fasta, e.message))

        self._sequence_rows = rows
        self._sequence = ''.join(rows)

    def _parse_sequence(self):
        """
        This method is used to parse the sequence and hash it on sequence identifier.
        It also stores the length of the longest sequence id
        """
        seq_hash = dict()
        seq_id = None

        # Python adaptation: Pre-compile regex
        re_1 = re.compile(r'^\s*$')
        re_2 = re.compile(r'^>(\S+)')

        for row in self._sequence_rows:

            # Ignore blank lines
            if re_1.search(row):
                continue

            res_re_2 = re_2.search(row)
            if res_re_2:
                seq_id = res_re_2.group(1)

                if seq_id in seq_hash:
                    sys.exit('FATAL: Sequence identifiers must be unique. Your fasta file contains two sequences with the same id (%s)' % seq_id)

                # Store the max length of seq name, use this later when printing in ascii
                if self._max_seqname is None or len(seq_id) > self._max_seqname:
                    self._max_seqname = len(seq_id)

            else:
                if not seq_id:
                    sys.exit("FATAL: Unrecognised format of fasta file. Each sequence must have a header line in the format '>identifier  <optional description>'")
                row = row.rstrip()
                seq_hash[seq_id] = row

        self._seq_hash = seq_hash

    def search(self):
        """
        The main method on the object. Performs a C<hmmscan> search using the supplied sequence and the specified HMM library.
        """

        if not self._sequence:
            sys.exit('FATAL: no sequence given; set the search parameters before calling "search"')

        AllResults = dict()
        pfamB = None
        firstResult = None

        # Python adaptation: Pre-compile regex.
        re_1 = re.compile(r'Pfam\-B')

        for hmmlib in self._hmmlib:
            hmmscan_cut_off = list()
            seq_evalue = None
            dom_evalue = None
            if not re_1.search(hmmlib):
                hmmscan_cut_off = self._hmmscan_cutoff
            else:
                pfamB = True
                seq_evalue = 0.001
                dom_evalue = 0.001

                # It's a pfamB search so use some default cut off values
                hmmscan_cut_off.append('-E %s' % seq_evalue)
                hmmscan_cut_off.append('--domE %s' % dom_evalue)

            if self._cpu:
                params = ['hmmsearch', '--notextw', '--cpu', str(self._cpu),
                          ' '.join(hmmscan_cut_off),
                          os.path.join(self._dir, hmmlib), self._fasta]
            else:
                params = ['hmmsearch', '--notextw', ' '.join(hmmscan_cut_off),
                          os.path.join(self._dir, hmmlib), self._fasta]

            proc = subprocess.Popen(params, stdout=subprocess.PIPE,
                                    stderr=subprocess.STDOUT, encoding='utf-8')
            proc_out, proc_err = proc.communicate()
            proc_out_ascii = proc_out

            if proc_err:
                sys.exit('An error was encountered while running hmmsearch: %s' % proc_err)

            self._hmmresultIO = HMMResultsIO()
            self._all_results = self._hmmresultIO.parseMultiHMMER3(proc_out_ascii)
            self._all_results = self._convert_results_search_to_scan(self._all_results)

            if not re_1.search(hmmlib):

                if self._clan_overlap:
                    print('Resolve clan overlaps: off')
                else:
                    self._resolve_clan_overlap()

            # Determine which hits are significant
            for result in self._all_results:

                for unit in sorted(result.units, key=lambda x: x.seqFrom, reverse=False):
                    if not pfamB:
                        unit.sig = False
                        if result.seqs[unit.name].bits >= self._seqGA[unit.name]:
                            if unit.bits >= self._domGA[unit.name]:
                                unit.sig = True
            if firstResult:
                AllResults[self._all_results] = self._all_results
            else:
                firstResult = self._all_results

        # If more than one search, merge results into one object
        if len(AllResults) > 0:
            sys.exit('FATAL: multiple search functionality has not been implemented in Python')

        # print('<seq id> <alignment start> <alignment end> <envelope start> <envelope end> <hmm acc> <hmm name> <type> <hmm start> <hmm end> <hmm length> <bit score> <E-value> <significance> <clan>')

    def write_results(self, out, e_seq, e_dom, b_seq, b_dom):
        """
        Writes the results of the C<hmmscan> search. Takes a single argument, which can
        be an open filehandle or a filename. A fatal error is generated if a file of the
        given name already exists.
        """
        # if os.path.isfile(out):
        #     sys.exit('FATAL: output file "%s" already exists' % out)

        try:
            with open(out, 'w') as fh:

                # Python adaptation: Write the header.
                fh.write('#' * 84 + '\n')
                fh.write('# %s #\n' % "'||'''|,          '||'''|,  .|';                   ".center(80, ' '))
                fh.write('# %s #\n' % " ||   ||           ||   ||  ||                     ".center(80, ' '))
                fh.write('# %s #\n' % " ||...|' '||  ||`  ||...|' '||'   '''|.  '||),,(|, ".center(80, ' '))
                fh.write('# %s #\n' % " ||       `|..||   ||       ||   .|''||   || || || ".center(80, ' '))
                fh.write('# %s #\n' % ".||           ||  .||      .||.  `|..||. .||    ||.".center(80, ' '))
                fh.write('# %s #\n' % "           ,  |'                                   ".center(80, ' '))
                fh.write('# %s #\n' % "            ''                                     ".center(80, ' '))
                fh.write('# %s #\n' % ''.center(80, ' '))
                fh.write('# %s #\n' % 'This program is free software: you can redistribute it and/or modify'.center(80, ' '))
                fh.write('# %s #\n' % 'it under the terms of the GNU General Public License as published by'.center(80, ' '))
                fh.write('# %s #\n' % 'the Free Software Foundation, either version 3 of the License, or'.center(80, ' '))
                fh.write('# %s #\n' % '(at your option) any later version.'.center(80, ' '))
                fh.write('# %s #\n' % 'This program is distributed in the hope that it will be useful,'.center(80, ' '))
                fh.write('# %s #\n' % 'but WITHOUT ANY WARRANTY; without even the implied warranty of'.center(80, ' '))
                fh.write('# %s #\n' % 'MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the'.center(80, ' '))
                fh.write('# %s #\n' % 'GNU General Public License for more details.'.center(80, ' '))
                fh.write('# %s #\n' % 'You should have received a copy of the GNU General Public License'.center(80, ' '))
                fh.write('# %s #\n' % 'along with this program. If not, see <http://www.gnu.org/licenses/>.'.center(80, ' '))
                fh.write('# %s #\n' % ''.center(80, ' '))
                fh.write('# %s #\n' % 'PyPfam is a Python (lite) port of the PfamSearch scripts written in Perl.'.center(80, ' '))
                fh.write('# %s #\n' % 'Perl authors: Jaina Mistry (jm14@sanger.ac.uk), John Tate (jt6@sanger.ac.uk), '.center(80, ' '))
                fh.write('# %s #\n' % 'Rob Finn (rdf@sanger.ac.uk)'.center(80, ' '))
                fh.write('# %s #\n' % 'Python authors: Aaron Mussig (a.mussig@uq.edu.au)'.center(80, ' '))
                fh.write('# %s #\n' % ''.center(80, '='))
                fh.write('# %-25s %s\n' % ('Execution time:', datetime.now()))
                fh.write('# %-25s %s\n' % ('Query sequence file:', self._fasta))
                fh.write('# %-25s %s\n' % ('CPU threads:', self._cpu))
                fh.write('# %-25s %s\n' % ('Searching against:', os.path.join(self._dir, 'Pfam-A.hmm')))
                fh.write('# %-25s %s\n' % ('Resolve clan overlaps:', 'on'))
                fh.write('# %s #\n' % ''.center(80, '='))
                fh.write('# <seq id> <alignment start> <alignment end> <envelope start> <envelope end> <hmm acc> <hmm name> <type> <hmm start> <hmm end> <hmm length> <bit score> <E-value> <significance> <clan> \n')

                for result in self._all_results:
                    self._hmmresultIO.write_ascii_out(result, fh, self, e_seq, e_dom, b_seq, b_dom)

        except IOError as e:
            sys.exit('FATAL: Can\'t write to your output file "%s": %s' % (out, e.message))

    def _resolve_clan_overlap(self):
        """
        Resolves overlaps between clans.
        """

        no_clan_overlap = list()
        for result in self._all_results:
            new = result.remove_overlaps_by_clan(self._clanmap, self._nested)
            no_clan_overlap.append(new)

        self._all_results = no_clan_overlap
