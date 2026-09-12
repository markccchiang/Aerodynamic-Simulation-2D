"""Sphinx configuration for the aerosim theory documentation."""

from __future__ import annotations

import os
import sys

# Put src/ on sys.path so the matplotlib plot_directive blocks below can
# ``from aerosim import ...`` and call the live solver to generate figures.
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

project = "aerosim"
author = "Cheng-Chin Chiang"
copyright = "%Y, Cheng-Chin Chiang"

extensions = [
    "sphinx.ext.mathjax",
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
    "sphinx_rtd_theme",
    "matplotlib.sphinxext.plot_directive",
]

# matplotlib plot_directive: PNG only, no source link, no format chooser —
# the figures should read as illustrations, not interactive widgets.
plot_include_source = False
plot_html_show_source_link = False
plot_html_show_formats = False
plot_formats = [("png", 120)]
plot_rcparams = {
    "figure.dpi": 120,
    "savefig.dpi": 120,
    "savefig.bbox": "tight",
    "font.size": 9,
    "axes.titlesize": 10,
    "axes.labelsize": 9,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "figure.facecolor": "white",
}
plot_apply_rcparams = True

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]
html_title = "aerosim — Aerodynamic Theory"
html_theme_options = {
    "navigation_depth": 3,
    "collapse_navigation": False,
    "sticky_navigation": True,
    "titles_only": False,
}

# Make cross-references to NumPy / SciPy resolve in docstrings.
intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "scipy": ("https://docs.scipy.org/doc/scipy/", None),
}

# Render math with MathJax (loaded from a CDN at build time).
mathjax3_config = {
    "tex": {
        "inlineMath": [["\\(", "\\)"], ["$", "$"]],
        "displayMath": [["\\[", "\\]"], ["$$", "$$"]],
    },
}

# So that source-block headings line up with the panel-method conventions used
# throughout the rest of the codebase.
rst_prolog = """
.. |Vinf| replace:: :math:`V_\\infty`
.. |alpha| replace:: :math:`\\alpha`
.. |sigma| replace:: :math:`\\sigma`
.. |gamma| replace:: :math:`\\gamma`
.. |Cp| replace:: :math:`C_p`
.. |Cl| replace:: :math:`C_\\ell`
.. |Cd| replace:: :math:`C_d`
.. |Cm| replace:: :math:`C_m`
.. |Re| replace:: :math:`\\mathrm{Re}`
"""
