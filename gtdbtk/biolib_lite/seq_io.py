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


import gzip
import os
import sys
import traceback

from .exceptions import InputFileError

protein_bases = {'a', 'r', 'n', 'd', 'c', 'q', 'e', 'g', 'h',
                 'i', 'l', 'k', 'm', 'f', 'p', 's', 't', 'w', 'y', 'v'}
nucleotide_bases = {'a', 'c', 'g', 't'}
insertion_bases = {'-', '.'}


def read_fasta(fasta_file, keep_annotation=False):
    """Read sequences from fasta file.

    Parameters
    ----------
    fasta_file : str
        Name of fasta file to read.
    keep_annotation : boolean
        Determine is sequence id should contain annotation.

    Returns
    -------
    dict : dict[seq_id] -> seq
        Sequences indexed by sequence id.
    """

    if not os.path.exists(fasta_file):
        raise InputFileError('Input file %s does not exist.' % fasta_file)

    if os.stat(fasta_file).st_size == 0:
        return {}

    try:

        if fasta_file.endswith('.gz'):
            file_f, file_mode = gzip.open, 'rt'
        else:
            file_f, file_mode = open, 'r'

        seqs = {}
        with file_f(fasta_file, file_mode) as f:

            for line in f.readlines():
                # skip blank lines
                if not line.strip():
                    continue

                if line[0] == '>':
                    if keep_annotation:
                        seq_id = line[1:-1]
                    else:
                        seq_id = line[1:].split(None, 1)[0]

                    seqs[seq_id] = []
                else:
                    seqs[seq_id].append(line.strip())

        for seq_id, seq in seqs.items():
            seqs[seq_id] = ''.join(seq).replace(' ', '')
    except:
        print(traceback.format_exc())
        print()
        print("[Error] Failed to process sequence file: " + fasta_file)
        sys.exit(1)

    return seqs


def read_fasta_seq(fasta_file, keep_annotation=False):
    """Generator function to read sequences from fasta file.

    This function is intended to be used as a generator
    in order to avoid having to have large sequence files
    in memory. Input file may be gzipped.

    Example:
    seq_io = SeqIO()
    for seq_id, seq in seq_io.read_fasta_seq(fasta_file):
        print seq_id
        print seq

    Parameters
    ----------
    fasta_file : str
        Name of fasta file to read.
    keep_annotation : boolean
        Determine if annotation string should be returned.

    Yields
    ------
    list : [seq_id, seq, [annotation]]
        Unique id of the sequence followed by the sequence itself,
        and the annotation if keep_annotation is True.
    """

    if not os.path.exists(fasta_file):
        raise InputFileError('Input file %s does not exist.' % fasta_file)

    if os.stat(fasta_file).st_size == 0:
        pass

    try:
        open_file = open
        mode = 'r'
        if fasta_file.endswith('.gz'):
            open_file = gzip.open
            mode = 'rb'

        seq_id = None
        annotation = None
        seq = None
        with open_file(fasta_file, mode) as f:

            for line in f.readlines():
                if isinstance(line, bytes):
                    line = line.decode()
                # skip blank lines
                if not line.strip():
                    continue

                if line[0] == '>':
                    if seq_id is not None:
                        if keep_annotation:
                            yield seq_id, ''.join(seq).replace(' ', ''), annotation
                        else:
                            yield seq_id, ''.join(seq).replace(' ', '')

                    line_split = line[1:-1].split(None, 1)
                    if len(line_split) == 2:
                        seq_id, annotation = line_split
                    else:
                        seq_id = line_split[0]
                        annotation = ''
                    seq = []
                else:
                    seq.append(line.strip())

        # report last sequence
        if seq:
            if keep_annotation:
                yield seq_id, ''.join(seq).replace(' ', ''), annotation
            else:
                yield seq_id, ''.join(seq).replace(' ', '')
    except GeneratorExit:
        pass
    except:
        print(traceback.format_exc())
        print()
        print("[Error] Failed to process sequence file: " + fasta_file)
        sys.exit(1)


def read_seq(seq_file, keep_annotation=False):
    """Generator function to read sequences from fasta/q file.

    This function is intended to be used as a generator
    in order to avoid having to have large sequence files
    in memory. Input file may be gzipped and in either
    fasta or fastq format. It is slightly more efficient
    to directly call read_fasta_seq() or read_fastq_seq()
    if the type of input file in known.

    Example:
    seq_io = SeqIO()
    for seq_id, seq in seq_io.read_seq(fasta_file):
        print seq_id
        print seq

    Parameters
    ----------
    seq_file : str
        Name of fasta/q file to read.
    keep_annotation : boolean
        Determine if annotation string should be returned.

    Yields
    ------
    list : [seq_id, seq, [annotation]]
        Unique id of the sequence followed by the sequence itself,
        and the annotation if keep_annotation is True.
    """

    if seq_file.endswith(('.fq.gz', '.fastq.gz', '.fq', '.fq.gz')):
        raise Exception("Cannot read FASTQ files.")
        # for rtn in read_fastq_seq(seq_file):
        #     yield rtn
    else:
        for rtn in read_fasta_seq(seq_file, keep_annotation):
            yield rtn


def write_fasta(seqs, fasta_file, wrap=80):
    """Write sequences to a fasta file.

    Parameters
    ----------
    seqs : dict[seq_id] -> seq
        Sequences indexed by sequence id.
    fasta_file : str
        Path to write the sequences to.
    wrap: int
        Number of AA/NT before the line is wrapped.
    """
    with open(fasta_file, 'w') as f:
        for gid, gseq in seqs.items():
            f.write('>{}\n'.format(gid))
            for i in range(0, len(gseq), wrap):
                f.write('{}\n'.format(gseq[i:i + wrap]))
