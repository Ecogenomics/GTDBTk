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

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2015'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'


class BioLibError(Exception):
    def __init__(self, message):
        Exception.__init__(self, message)


class BioLibFileNotFound(BioLibError):
    """ Raised when a file is not found. """

    def __init__(self, message=''):
        super(BioLibFileNotFound, self).__init__(message)


class BioLibDirNotFound(BioLibError):
    """ Raised when a directory is not found. """

    def __init__(self, message=''):
        super(BioLibDirNotFound, self).__init__(message)


class BioLibIOException(BioLibError):
    """ Raised when an IO exception is encountered. """

    def __init__(self, message=''):
        super(BioLibIOException, self).__init__(message)


class InputFileError(BioLibError):
    """Raised when an input file error is encountered."""

    def __init__(self, message=''):
        super(InputFileError, self).__init__(message)
