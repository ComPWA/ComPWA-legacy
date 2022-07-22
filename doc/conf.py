# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

import os

# -- Project information -----------------------------------------------------

project = "ComPWA"
copyright = "2020, The ComPWA Team"
author = "The ComPWA Team"


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "breathe",
    "exhale",
    "sphinx.ext.githubpages",
    "sphinx.ext.mathjax",
    "sphinx.ext.intersphinx",
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = [
    "_build",
    "Thumbs.db",
    ".DS_Store",
    "api/program_listing_file_doc_main.md.rst",
]

# -- Extension settings ------------------------------------------------------

# Setup Sphinx
primary_domain = "cpp"
highlight_language = "cpp"
html_logo = "images/logo-small.png"
html_theme = "sphinx_book_theme"
html_theme_options = {
    "repository_url": "https://github.com/ComPWA/ComPWA",
    "use_repository_button": True,
}

# Setup Breathe
breathe_projects = {"ComPWA": "./doxyoutput/xml"}
breathe_default_project = "ComPWA"

# Setup Exhale
exhale_args = {
    # These arguments are required
    "containmentFolder": "./api",
    "rootFileName": "index.rst",
    "rootFileTitle": "API",
    "doxygenStripFromPath": "..",
    # Suggested optional arguments
    "createTreeView": True,
    "exhaleExecutesDoxygen": True,
    "exhaleUseDoxyfile": True,
}

# Setup intersphinx
intersphinx_mapping = {
    "pycompwa": ("https://compwa.github.io/", None),
}
