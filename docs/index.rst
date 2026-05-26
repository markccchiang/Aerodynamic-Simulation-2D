aerosim — Aerodynamic Theory
============================

`aerosim` is a 2D airfoil aerodynamics simulator built from scratch in
Python: a **Hess–Smith panel method** for the inviscid flow plus an
**integral boundary-layer correction** for the viscous profile drag, with
an optional **two-way viscous–inviscid coupling**.

This document explains the aerodynamic theory the simulator implements:
the governing equations, why a panel method is the natural way to solve
them, what the boundary-layer correction adds, and how the two are
coupled.  The chapters below are largely self-contained but build up in
the order the solver itself executes.

.. toctree::
   :maxdepth: 2
   :caption: Theory

   theory/overview
   theory/governing_equations
   theory/panel_method
   theory/boundary_layer
   theory/coupling
   theory/validation

.. toctree::
   :maxdepth: 1
   :caption: Reference

   theory/symbols
   theory/references


Building these docs
-------------------

.. code-block:: bash

   uv add --dev sphinx sphinx_rtd_theme    # one-time
   cd docs && uv run make html             # output in docs/_build/html/index.html


Indices
-------

* :ref:`genindex`
* :ref:`search`
