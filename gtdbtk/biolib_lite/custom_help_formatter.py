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
__credits__ = ['Donovan Parks', 'Connor Skennerton']
__license__ = 'GPL3'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'

import argparse
import os
import tempfile


class ChangeTempAction(argparse.Action):
    """Action for changing the directory used for temporary files.

    Example:
        <parse>.add_argument('--tmpdir', action=ChangeTempAction, default=tempfile.gettempdir(), help="specify alternative directory for temporary files")
    """

    def __call__(self, parser, namespace, values, option_string=None):
        if os.path.isdir(values):
            tempfile.tempdir = values
            setattr(namespace, self.dest, values)
        else:
            raise argparse.ArgumentTypeError(
                'The value of %s must be a valid directory' % option_string)


class CustomHelpFormatter(argparse.HelpFormatter):
    """Provide a customized format for help output.

    http://stackoverflow.com/questions/9642692/argparse-help-without-duplicate-allcaps
    """

    def _get_help_string(self, action):
        """Place default value in help string."""
        h = action.help

        # Remove any formatting used for Sphinx argparse hints.
        h = h.replace('``', '')

        if '%(default)' not in action.help:
            if action.default != '' and action.default != [] and \
                    action.default is not None and \
                    not isinstance(type(action.default), bool):
                if action.default is not argparse.SUPPRESS:
                    defaulting_nargs = [
                        argparse.OPTIONAL, argparse.ZERO_OR_MORE]

                    if action.option_strings or action.nargs in defaulting_nargs:
                        if '\n' in h:
                            lines = h.splitlines()
                            lines[0] += ' (default: %(default)s)'
                            h = '\n'.join(lines)
                        else:
                            h += ' (default: %(default)s)'
            return h

    def _format_action_invocation(self, action):
        """Removes duplicate ALLCAPS with positional arguments."""
        if not action.option_strings:
            default = self._get_default_metavar_for_positional(action)
            metavar, = self._metavar_formatter(action, default)(1)
            return metavar

        else:
            parts = []

            # if the Optional doesn't take a value, format is:
            #    -s, --long
            if action.nargs == 0:
                parts.extend(action.option_strings)

            # if the Optional takes a value, format is:
            #    -s ARGS, --long ARGS
            else:
                default = self._get_default_metavar_for_optional(action)
                args_string = self._format_args(action, default)
                for option_string in action.option_strings:
                    parts.append(option_string)

                return '%s %s' % (', '.join(parts), args_string)

            return ', '.join(parts)

    def _get_default_metavar_for_optional(self, action):
        return action.dest.upper()

    def _get_default_metavar_for_positional(self, action):
        return action.dest
