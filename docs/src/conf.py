import os
import sys
from datetime import datetime
import re

sys.path.insert(0, os.path.abspath('../..'))
from gtdbtk import __author__, __title__, __maintainer__, __url__

# Configuration file for the Sphinx documentation builder.
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------

project = __title__
copyright = f'{datetime.now().year}, {__maintainer__}'
author = __author__

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'sphinx'

github_url = __url__

def get_version():
    with open(os.path.abspath('../../gtdbtk/__init__.py')) as f:
        for line in f:
            if line.startswith('__version__'):
                return re.search(r'["\'](.+?)["\']', line).group(1)
    return 'unknown'

release = get_version()
version = '.'.join(release.split('.')[:2])

# Make version info available for substitution in RST files
rst_epilog = f"""
.. |release| replace:: {release}
.. |version| replace:: {version}
"""

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ['sphinxarg.ext', 'sphinx.ext.napoleon', 'sphinx.ext.autodoc', 'linuxdoc.rstFlatTable',
              'recommonmark', 'sphinx_sitemap', 'nbsphinx','matplotlib.sphinxext.plot_directive']

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []

source_suffix = {
    '.rst': 'restructuredtext',
    '.txt': 'markdown',
    '.md': 'markdown',
}

html_extra_path = ['robots.txt']

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.

html_theme = 'sphinx_rtd_theme'

html_theme_options = {
    'canonical_url': '',
    # 'analytics_id': 'UA-84847737-2',  # Provided by Google in your dashboard
    'logo_only': True,
    'display_version': True,
    'prev_next_buttons_location': 'bottom',
    'style_external_links': True,
    # Toc options
    'collapse_navigation': True,
    'sticky_navigation': True,
    'navigation_depth': 4,
    'includehidden': True,
    'titles_only': False
}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

html_context = {
    'css_files': [
        '_static/theme_overrides.css',  # override wide tables in RTD theme
        ],
    'display_github': True,
    'github_user': 'Ecogenomics',
    'github_repo': 'GTDBTk',
    'github_version': 'master/docs/src/'
     }

html_js_files = [
    'custom.js'
]

html_logo = '_static/GTDBTk.svg'

# Sitemap settings
html_baseurl = 'https://ecogenomics.github.io/GTDBTk/'
sitemap_url_scheme = "{link}"


def get_version():
    with open(os.path.abspath('../../gtdbtk/__init__.py')) as f:
        for line in f:
            if line.startswith('__version__'):
                return re.search(r'["\'](.+?)["\']', line).group(1)
    return 'unknown'

release_parsed = get_version()
version_parsed = '.'.join(release.split('.')[:2])

# Make version info available for substitution in RST files
rst_epilog = f"""
.. |release| replace:: {release_parsed}
.. |version| replace:: {version_parsed}
"""
