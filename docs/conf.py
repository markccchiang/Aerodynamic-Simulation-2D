"""Sphinx configuration for the aerosim theory documentation."""

from __future__ import annotations

project = "aerosim"
author = "aerosim contributors"
copyright = "%Y, aerosim contributors"

extensions = [
    "sphinx.ext.mathjax",
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# Default theme; swap for "furo" or "sphinx_rtd_theme" if installed.
html_theme = "alabaster"
html_static_path = ["_static"]
html_title = "aerosim — Aerodynamic Theory"

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
