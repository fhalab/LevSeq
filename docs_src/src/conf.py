import os
import sys
from datetime import datetime

sys.path.insert(0, os.path.abspath('../..'))
__author__ = 'Yueming,Long'
__version__ = '0.0.1'
__title__ = 'LevSeq'
__url__ = 'https://github.com/fhalab/LevSeq'

# Configuration file for the Sphinx documentation builder.
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
import sphinx_rtd_theme

html_theme = "sphinx_rtd_theme"
html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]

project = __title__
copyright = f'{datetime.now().year}, {__author__}'
author = __author__

# The full version, including alpha/beta/rc tags
release = __version__
version = __version__

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'sphinx'

github_url = __url__

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ['sphinx.ext.napoleon', 'sphinx.ext.autodoc',
               'recommonmark', 'sphinx_sitemap', 'nbsphinx',
              'sphinxcontrib.jquery']

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []

source_suffix = {
    '.rst': 'restructuredtext',
    '.txt': 'markdown',
    '.md': 'markdown'
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
    # 'css_files': [
    #     '_static/theme_overrides.css',  # override wide tables in RTD theme
    #     ],
    'display_github': True,
    'github_user': 'fhalab',
    'github_repo': 'LevSeq',
    'github_version': 'main/docs_src/src/'
     }

html_js_files = [
    'custom.js'
]

html_logo = '_static/logo.png'

# Sitemap settings
html_baseurl = 'https://github.com/fhalab/LevSeq/'
sitemap_url_scheme = "{link}"

nbsphinx_allow_errors = True

# Favicon
html_favicon = '_static/logo.png'
