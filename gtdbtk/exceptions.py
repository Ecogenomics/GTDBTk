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
    pass


class GenomeNameInvalid(GTDBTkException):
    """ Thrown when a genome name contains characters which are not supported. """
    pass


class GenomeBatchfileMalformed(GTDBTkException):
    """ Thrown when the format of the genome batchfile is malformed. """
    pass


class NoGenomesFound(GTDBTkException):
    """ Thrown when no input genomes are found in a directory. """
    pass


class ReferenceFileMalformed(GTDBTkException):
    """ Thrown when a reference file is malformed. """
    pass


class GenomeMarkerSetUnknown(GTDBTkException):
    """ Thrown when the genome marker set is unknown (i.e. not ar122, or bac120). """
    pass


class FileNotFound(GTDBTkException):
    """ Thrown when a file is not found. """
    pass


class DirNotFound(GTDBTkException):
    """ Thrown when a directory is not found. """


class ProdigalException(GTDBTkException):
    """ Thrown when Prodigal returns a non-zero exit code. """
