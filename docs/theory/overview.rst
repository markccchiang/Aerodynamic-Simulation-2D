Overview
========

What the simulator computes
---------------------------

For a given 2D airfoil shape and angle of attack |alpha|, `aerosim`
predicts the steady, incompressible flow around the section and returns:

* the **surface pressure distribution** :math:`C_p(s)` along the chord;
* the **lift, drag, and pitching-moment coefficients** :math:`C_\ell`,
  :math:`C_d`, :math:`C_m`;
* an off-body **velocity field** suitable for streamline plotting;
* and, optionally, a **viscous boundary-layer profile** with the local
  momentum thickness :math:`\theta(s)`, displacement thickness
  :math:`\delta^*(s)`, shape factor :math:`H(s)`, transition location,
  and a profile-drag estimate.

The model is two-dimensional and assumes **incompressible, irrotational,
attached** (or only mildly separated) flow at small angles of attack —
the classical operating regime for an airfoil section well below stall.


The solve pipeline
------------------

The execution order mirrors the physics:

.. math::

    \text{Geometry}
    \;\xrightarrow{\text{Hess--Smith}}\;
    \underbrace{(C_p,\,V_t,\,C_\ell,\,C_m)}_{\text{inviscid}}
    \;\xrightarrow{\text{integral BL}}\;
    \underbrace{(\theta,\,\delta^*,\,H,\,x_\text{tr},\,C_{d,\text{visc}})}_{\text{viscous correction}}

When :ref:`two-way coupling <coupling>` is enabled, the boundary-layer
displacement is fed back into the panel solve as a wall-transpiration
boundary condition and the loop is iterated to convergence, so that
:math:`C_\ell` and :math:`C_p` themselves react to viscosity (the
**viscous decambering** effect).

In picture form:

.. code-block:: text

       ┌──────────────────┐
       │     Geometry     │   nodes, panels, cosine spacing
       └────────┬─────────┘
                │
                ▼
       ┌──────────────────┐
       │   Hess–Smith     │   (N+1) linear system
       │  panel solver    │   sources σⱼ  +  shared vortex γ
       └────────┬─────────┘
                │   Cp, Vt, Cl, Cm
                ▼
       ┌──────────────────┐                      ┌───────────────────────┐
       │  Boundary layer  │  edge velocity Ue    │   Wall-transpiration  │
       │  Thwaites→Michel │────────────────────▶ │   v_n = d(Ue·δ*)/ds   │
       │  →Head + LT      │                      └───────────┬───────────┘
       └────────┬─────────┘                                  │
                │   θ, δ*, H, x_tr, Cd_visc                  │  (couple=True only:
                │                                            │   feed back into the
                ▼                                            │   flow-tangency RHS
       ┌──────────────────┐                                  │   and re-solve)
       │     Solution     │ ◀────────────────────────────────┘
       └──────────────────┘

The dashed feedback path is what ``couple=True`` switches on; without
it, the boundary-layer block is a pure one-way post-process.


Why two solvers?
----------------

The inviscid panel solver alone is fast and accurate for **lift** and
the **pressure distribution** at small angles, but it predicts
**zero drag** (a result known as **d'Alembert's paradox**, derived in
:doc:`governing_equations`). It also has no notion of stall.

Real drag is dominated by what happens in the thin boundary layer at
the surface:

* **viscous skin friction** along the wetted surface;
* a small **form drag** due to the boundary-layer thickening the
  effective shape;
* and — once separation grows — pressure drag from a wake that the
  inviscid model cannot represent.

A full Navier–Stokes solver would resolve all of this, at high cost. An
**integral boundary-layer method** is a much cheaper compromise: it
reduces the boundary-layer PDE to one or two ODEs in arc length that
ride along the inviscid surface-velocity distribution. It gets profile
drag and transition right to a few percent for attached flow, and tells
you honestly when it can no longer cope (separation, post-stall).

The chapters that follow develop the two halves of the model and the
optional coupling between them:

* :doc:`governing_equations` — incompressible potential flow, the
  Kutta condition, d'Alembert.
* :doc:`panel_method` — discretising the surface, source-plus-vortex
  singularities, the :math:`(N+1)` linear system, integration to
  :math:`C_\ell` / :math:`C_m`.
* :doc:`boundary_layer` — Thwaites, Michel transition, Head's
  entrainment method, Squire–Young.
* :doc:`coupling` — wall-transpiration model, under-relaxed direct
  iteration, where it stops working.
* :doc:`validation` — what the test suite checks and against what.
