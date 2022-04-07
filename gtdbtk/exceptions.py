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


class GTDBTkException(Exception):
    """ Base exception for all GTDB-Tk exceptions thrown in this project. """

    def __init__(self, message=''):
        Exception.__init__(self, message)


class GTDBTkExit(Exception):
    """Raised when GTDB-Tk is to quietly exit."""

    def __init__(self, message=''):
        Exception.__init__(self, message)


class GenomeNameInvalid(GTDBTkException):
    """ Thrown when a genome name contains characters which are not supported. """

    def __init__(self, message=''):
        GTDBTkException.__init__(self, message)


class GenomeBatchfileMalformed(GTDBTkException):
    """ Thrown when the format of the genome batchfile is malformed. """

    def __init__(self, message=''):
        GTDBTkException.__init__(self, message)


class NoGenomesFound(GTDBTkException):
    """ Thrown when no input genomes are found in a directory. """

    def __init__(self, message=''):
        GTDBTkException.__init__(self, message)


class ReferenceFileMalformed(GTDBTkException):
    """ Thrown when a reference file is malformed. """

    def __init__(self, message=''):
        GTDBTkException.__init__(self, message)


class GenomeMarkerSetUnknown(GTDBTkException):
    """ Thrown when the genome marker set is unknown (i.e. not ar53, or bac120). """

    def __init__(self, message=''):
        GTDBTkException.__init__(self, message)


class InconsistentGenomeBatch(GTDBTkException):
    """ Thrown when number of genomes in the identify directory is different than the number of genomes to process. """

    def __init__(self, message=''):
        GTDBTkException.__init__(self, message)


class FileNotFound(GTDBTkException):
    """ Thrown when a file is not found. """

    def __init__(self, message=''):
        GTDBTkException.__init__(self, message)


class DirNotFound(GTDBTkException):
    """ Thrown when a directory is not found. """

    def __init__(self, message=''):
        GTDBTkException.__init__(self, message)


class ProdigalException(GTDBTkException):
    """ Thrown when Prodigal returns a non-zero exit code. """

    def __init__(self, message=''):
        GTDBTkException.__init__(self, message)


class MSAMarkerLengthMismatch(GTDBTkException):
    """ Thrown when an MSA length does not equal the length of the markers. """

    def __init__(self, message=''):
        GTDBTkException.__init__(self, message)


class MSAMaskLengthMismatch(GTDBTkException):
    """ Thrown when an MSA length does not equal the length of the mask. """

    def __init__(self, message=''):
        GTDBTkException.__init__(self, message)


class PplacerException(GTDBTkException):
    """ Thrown whenever an error is encountered while running pplacer. """

    def __init__(self, message=''):
        GTDBTkException.__init__(self, message)


class TogException(GTDBTkException):
    """ Thrown whenever an error is encountered while running tog. """

    def __init__(self, message=''):
        GTDBTkException.__init__(self, message)


class FastANIException(GTDBTkException):
    """ Thrown when an exception is encountered while running FastANI. """

    def __init__(self, message=''):
        GTDBTkException.__init__(self, message)


class FastTreeException(GTDBTkException):
    """ Thrown when an exception is encountered while running FastTree. """

    def __init__(self, message=''):
        GTDBTkException.__init__(self, message)


class GTDBTkTestFailure(GTDBTkException):
    """ Thrown when the GTDBTk user test suite fails. """

    def __init__(self, message=''):
        GTDBTkException.__init__(self, message)


class GTDBTkDataPathUndefined(GTDBTkException):
    """ Thrown when the GTDBTK_DATA_PATH environment variable is undefined. """

    def __init__(self, message=''):
        GTDBTkException.__init__(self, message)

class GTDBTkArgsParsingConflict(GTDBTkException):
    """ Thrown when the arguments are conflicting and or missing.  """

    def __init__(self, message=''):
        GTDBTkException.__init__(self, message)
