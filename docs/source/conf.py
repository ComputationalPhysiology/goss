# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html
# -- Path setup --------------------------------------------------------------
# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# -- Project information -----------------------------------------------------
import goss

project = "goss"
copyright = "2022, Simula Research Laboratory"
author = "Simula Research Laboratory"

# The full version, including alpha/beta/rc tags


release = goss.__version__


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.mathjax",
    "sphinx.ext.viewcode",
    "sphinx.ext.napoleon",
    "sphinx.ext.githubpages",
    "sphinxcontrib.bibtex",
    "sphinx_math_dollar",
    "myst_parser",
]

bibtex_bibfiles = ["refs.bib"]

html_theme_options = {
    "toc_title": "goss",
    "repository_url": "https://github.com/ComputationalPhysiology/goss",
    "use_issues_button": True,
    "use_repository_button": True,
    "extra_navbar": '<a href="https://www.simula.no/research/projects/department-computational-physiology">Computational Physiology at Simula</a>',
}


html_logo = "logo.png"
html_favicon = "favicon.ico"

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
# exclude_patterns = []


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "sphinx_book_theme"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]

# Make sure we can reference figures with numbers
numfig = True
