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
__copyright__ = 'Copyright 2014'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'

import logging
import ntpath
import os

from .exceptions import BioLibFileNotFound, BioLibDirNotFound, BioLibIOException


def is_float(s):
    """Check if a string can be converted to a float.

    Parameters
    ----------
    s : str
        String to evaluate.

    Returns
    -------
    boolean
        True if string can be converted, else False.
    """
    try:
        float(s)
    except ValueError:
        return False

    return True


def check_file_exists(input_file):
    """Assert that a file exists.

    Parameters
    ----------
    input_file : str
        The path to file being checked.

    Returns
    -------
    bool
        True if the directory exists.

    Raises
    ------
    BioLibFileNotFound
        If the file doesn't exist.
    """
    if not os.path.exists(input_file) or not os.path.isfile(input_file):
        logger = logging.getLogger('timestamp')
        logger.error('Input file does not exist: ' + input_file + '\n')
        raise BioLibFileNotFound('Input file does not exist: ' + input_file + '\n')
    return True


def check_dir_exists(input_dir):
    """Assert that a directory exists.

    Parameters
    ----------
    input_dir : str
        The path to the directory being checked.

    Returns
    -------
    bool
        True if the directory exists.

    Raises
    ------
    BioLibDirNotFound
        If the directory doesn't exist.
    """
    if not os.path.exists(input_dir) or not os.path.isdir(input_dir):
        logger = logging.getLogger('timestamp')
        logger.error('Input directory does not exist: ' + input_dir + '\n')
        raise BioLibDirNotFound('Input directory does not exist: ' + input_dir + '\n')
    return True


def make_sure_path_exists(path):
    """Create directory if it does not exist.

    Parameters
    ----------
    path : str
        The path to the directory which should be created.

    Returns
    -------
    bool
        True if the path exists.

    Raises
    ------
    BioLibIOException
        If an error was encountered while creating the directory.
    """
    if not path:
        # lack of a path qualifier is acceptable as this
        # simply specifies the current directory
        return True
    elif os.path.isdir(path):
        return True

    try:
        os.makedirs(path)
        return True
    except OSError:
        logger = logging.getLogger('timestamp')
        logger.error('Specified path could not be created: ' + path + '\n')
        raise BioLibIOException('Specified path could not be created: ' + path + '\n')


def remove_extension(filename, extension=None):
    """Remove extension from filename. A specific extension can be specified,
    otherwise the extension is taken as all characters after the last period.

    Parameters
    ----------
    filename : str
        The name of the file in which the extension will be removed.
    extension : str, optional
        The extension to be removed.

    Returns
    -------
    str
        The name of the file with the extension removed.
    """
    f = ntpath.basename(filename)

    if extension and f.endswith(extension):
        f = f[0:f.rfind(extension)]
    else:
        f = os.path.splitext(f)[0]

    if f[-1] == '.':
        f = f[0:-1]

    return f
