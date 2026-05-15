import os
import sys
sys.path.insert(0, os.path.abspath('..'))
sys.path.insert(0, os.path.abspath('../kdellipspy'))

project = 'KDEllipsPy'
copyright = '2026, alexvillarroel'
author = 'alexvillarroel'

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.mathjax',
    'sphinx.ext.viewcode',
    'sphinx.ext.napoleon',
    'nbsphinx',
]

autodoc_mock_imports = [
    'numpy',
    'matplotlib',
    'pyproj',
    'scipy',
    'obspy',
    'neighpy',
]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store', '**.ipynb_checkpoints']

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

nbsphinx_execute = 'never'
