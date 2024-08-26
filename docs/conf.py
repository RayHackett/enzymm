# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys

sys.path.insert(
    0, os.path.realpath(os.path.join(__file__, os.path.pardir, os.path.pardir))
)
from . import __version__

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "matcher"
copyright = "2024, Raymund Hackett"
author = "Raymund Hackett"
version = __version__

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.intersphinx",
]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "alabaster"
html_static_path = ["_static"]

# -- Options for InterSphinx -------------------------------------------------

default_role = "py:obj"

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", (None, "python-inv.txt")),
    "pyjess": ("https://pyjess.readthedocs.io/en/latest/", None),
}
