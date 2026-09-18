.. _coupling:

Two-way viscous–inviscid coupling
=================================

Motivation: viscous decambering
-------------------------------

The uncoupled boundary-layer correction of :doc:`boundary_layer` gives
us a profile drag, but the inviscid surface pressure is left untouched.
In reality, the boundary layer **does** alter the pressure field: the
outer flow no longer sees the bare airfoil, it sees the airfoil
*plus* its displacement thickness :math:`\delta^*(s)`. That is a
slightly thicker, slightly de-cambered shape, which carries a little
less lift and has a softer suction peak — the familiar **viscous
decambering** effect.

Capturing this requires *feeding the boundary layer back* into the
panel solve. There are several classical strategies:

* **Direct coupling** — replace the body's no-penetration condition by
  a *wall-transpiration* condition that mimics the displacement effect,
  then iterate.  Robust in attached flow and through mild separation,
  but breaks down past stall (the **Goldstein singularity**).
* **Inverse coupling** — prescribe the displacement thickness and let
  the panel solver compute the edge velocity that produces it. Stable
  in massive separation but requires a different inviscid driver.
* **Semi-inverse coupling** — switch modes adaptively (e.g.,
  Veldman's method) so the same scheme works pre- and post-stall.

`aerosim` implements **direct coupling** with under-relaxation. It is
deliberately scoped to attached and mildly separated flow and reports
honestly (``Solution.converged``) when it cannot cope.


The transpiration model
-----------------------

The key trick is due to Lighthill: a viscous boundary layer of
displacement thickness :math:`\delta^*` looks, to the outer flow, like
an inviscid flow over a body that is *thicker* by :math:`\delta^*`. To
first order, that effect can be enforced as a small **blowing velocity
through the original surface**, normal to the wall:

.. math::
   :label: vn

   v_n(s) \;=\; \frac{1}{\rho_e}\,\frac{d(\rho_e U_e \delta^*)}{ds}
            \;\xrightarrow{\text{incompressible}}\;
            \frac{d(U_e\, \delta^*)}{ds}.

This is a 1D continuity statement: the mass that the displacement
thickness diverts from the wall has to come from somewhere, and in the
transpiration model it is *injected through the wall* at the rate
:math:`v_n`.  Equation :eq:`vn` is what
:func:`aerosim.boundary_layer.transpiration_velocity` evaluates.

The inviscid flow-tangency condition then becomes

.. math::
   :label: tang_visc

   \mathbf{V}\cdot \hat{\mathbf{n}}_\text{ext} \;=\; v_n(s),

instead of zero. In the discrete panel system this means adding
:math:`v_n` (sign-corrected by the exterior-side flag) to the RHS of the
:math:`N` flow-tangency equations — see
:func:`aerosim.panel._solve_panels` with ``vn`` supplied.


The iteration
-------------

With :eq:`vn` and :eq:`tang_visc` the inviscid and viscous problems are
mutually coupled (each needs the other's output). The simplest scheme is
**direct fixed-point iteration**:

.. math::
   :label: fixed_point

   \begin{aligned}
   \text{step 0:}\quad & v_n^{(0)} = 0
     \;\Rightarrow\; \text{inviscid }(\sigma,\gamma,V_t)^{(0)}, \\
   \text{step }k:\quad
     & \theta, \delta^*, H \;\leftarrow\; \text{BL}(U_e^{(k-1)};\, \mathrm{Re}), \\
     & \tilde v_n           \;\leftarrow\; d(U_e^{(k-1)}\delta^*)/ds, \\
     & v_n^{(k)}            \;\leftarrow\;
            (1-\omega)\, v_n^{(k-1)} + \omega\, \tilde v_n, \\
     & (\sigma,\gamma,V_t)^{(k)} \;\leftarrow\; \text{panel}(v_n^{(k)}),
   \end{aligned}

iterated until :math:`|C_\ell^{(k)} - C_\ell^{(k-1)}| < \varepsilon`.
Here :math:`\omega` is the **under-relaxation** factor.

The default settings used in `aerosim` are

.. math::

   \omega = 0.25, \qquad \varepsilon = 2 \times 10^{-5}, \qquad
   k_\text{max} = 80,

which converges in roughly 10–40 iterations across the usable range. A
larger :math:`\omega` (0.4) tends to oscillate under heavier loading
(high :math:`\alpha`, low :math:`\mathrm{Re}`), which is why the default
is conservative.


Practical fixes for cusped trailing edges
-----------------------------------------

Two numerical artifacts had to be tamed for the iteration to converge:

1. **TE cusp spikes.** A closed (cusped) trailing edge makes the
   inviscid :math:`U_e` drop sharply over the final ~0.3 % of chord,
   which spikes :math:`d(U_e \delta^*)/ds` and hence :math:`v_n` at the
   TE. Left raw, :math:`v_n^\text{max}` reaches O(10), which diverges
   the iteration immediately.

2. **Discrete noise.** The numerical derivative
   :math:`d(U_e\delta^*)/ds` is noisy on a finite panel grid.

`aerosim` therefore applies a **1-2-1 smoothing pass** to the mass-defect
quantity :math:`U_e \delta^*` and then **tapers the transpiration
velocity linearly to zero from** :math:`x/c = 0.95` **to the trailing
edge**:

.. math::
   v_n(x) \;\leftarrow\; v_n(x)\,
        \min\!\left(1,\; \frac{1 - x/c}{1 - 0.95}\right).

This is the same constant ``_TE_TAPER_X`` that governs the "separation
ignored inside the last 10% of chord" rule in :doc:`boundary_layer` —
both are dealing with the same cusped-TE artifact.


What changes when coupling is on
--------------------------------

Switching from ``solve(geom, alpha, re=Re)`` to
``solve(geom, alpha, re=Re, couple=True)`` makes three things differ:

* **Lift drops slightly.** Typical numbers from the validation suite
  for a NACA 2412 at :math:`\alpha = 5°,\,\mathrm{Re}=10^6`:
  :math:`C_\ell` falls from about 0.86 to 0.77 — the viscous decambering.
* **The suction peak softens.** The minimum of :math:`C_p` on the upper
  surface becomes less negative.
* **A small physical "pressure drag" appears.** Because transpiration
  makes the body slightly *leaky*, the pressure-integrated
  ``Solution.cd`` is no longer ~0. This is **expected**, not a bug:
  it is the inviscid form-drag contribution induced by the displacement
  body. The total profile drag is still
  ``Solution.cd_visc`` from Squire–Young.

The two effects are visible directly on the surface pressure plot:

.. plot::
   :alt: Inviscid vs viscously-coupled Cp on a NACA 2412.
   :caption: Surface :math:`C_p` for a NACA 2412 at α = 5°,
       Re = 10\ :sup:`6`. The coupled curve has a softer suction peak
       and a slightly higher pressure on the upper surface aft of mid-chord
       — the viscous decambering effect.

   import matplotlib.pyplot as plt

   from aerosim import naca4, Geometry, solve

   geom = Geometry(*naca4("2412"))
   inv = solve(geom, alpha_deg=5.0, re=1e6)
   vis = solve(geom, alpha_deg=5.0, re=1e6, couple=True)

   fig, ax = plt.subplots(figsize=(6.4, 3.6))
   ax.plot(inv.geom.xc, inv.cp, lw=1.2, label=f"inviscid  ($C_\\ell$ = {inv.cl:.3f})")
   ax.plot(vis.geom.xc, vis.cp, lw=1.2, ls="--",
           label=f"coupled ($C_\\ell$ = {vis.cl:.3f}, n_iter = {vis.n_iter})")
   ax.axhline(0, color="0.5", lw=0.6)
   ax.invert_yaxis()
   ax.set_xlabel("x / c"); ax.set_ylabel(r"$C_p$")
   ax.set_title("NACA 2412 surface pressure: inviscid vs coupled")
   ax.legend(frameon=False, fontsize=9)
   plt.tight_layout()

A look inside the iteration of equation :eq:`fixed_point` shows what
"converges in 10–40 iterations" actually means: the lift coefficient
drops in the first couple of steps, then the under-relaxation walks it
into the converged value.

.. plot::
   :alt: Lift-coefficient convergence history of the coupling iteration.
   :caption: Convergence of the under-relaxed direct iteration. The
       inviscid lift on iteration 0 falls quickly as the boundary-layer
       displacement feeds back, then settles to within
       :math:`2\times 10^{-5}` per step (the convergence tolerance).

   import numpy as np
   import matplotlib.pyplot as plt

   from aerosim import naca4, Geometry
   from aerosim.panel import _solve_panels, _COUPLE_RELAX
   from aerosim.boundary_layer import boundary_layer, transpiration_velocity

   geom = Geometry(*naca4("2412"))
   alpha, re, vinf = 5.0, 1e6, 1.0
   relax = _COUPLE_RELAX

   sol = _solve_panels(geom, alpha, vinf, vn=None)
   vn = np.zeros(geom.n)
   cls = [sol.cl]
   for _ in range(30):
       bl = boundary_layer(sol, re)
       vn = (1 - relax) * vn + relax * transpiration_velocity(sol, bl)
       sol = _solve_panels(geom, alpha, vinf, vn=vn)
       cls.append(sol.cl)

   fig, ax = plt.subplots(figsize=(5.6, 3.2))
   ax.plot(range(len(cls)), cls, "o-", ms=4, lw=1.2)
   ax.axhline(cls[-1], color="0.6", ls="--", lw=0.7)
   ax.set_xlabel("iteration k")
   ax.set_ylabel(r"$C_\ell^{(k)}$")
   ax.set_title("Under-relaxed direct coupling (ω = 0.25)")
   plt.tight_layout()

The lift loss grows as :math:`\mathrm{Re}` falls (thicker boundary
layer), and at :math:`\alpha = 0` on a symmetric airfoil it stays at
exactly zero — the coupling does not break the up/down symmetry that
would produce spurious lift. Both are part of the validation suite.


Where direct coupling stops working
-----------------------------------

Direct coupling has a fundamental obstacle in fully separated flow: the
**Goldstein singularity** [Gol1948]_. As separation approaches, the
inviscid system tries to extrapolate the boundary-layer growth
infinitely, and the iteration ceases to be a contraction. Symptoms:

* :math:`v_n` grows unboundedly near the separation point;
* the iteration oscillates or diverges;
* convergence reports ``False`` and the result should not be trusted.

Pushing through the singularity requires **inverse or semi-inverse**
coupling: invert the inviscid operator so that the displacement
thickness is the *input* to the panel solve, not the output. That is
explicitly out of scope here, and the solver reports
``Solution.converged = False`` for cases where it cannot make the
direct iteration close. The "Possible next steps" section of the
project's development guide (``guide/development.md``) lists
semi-inverse coupling as the natural next addition.


Where it fits in the code
-------------------------

The relevant entry points are:

* :func:`aerosim.panel.solve` — dispatches to ``_solve_coupled`` when
  ``couple=True`` (requires ``re=``).
* :func:`aerosim.panel._solve_panels` — single panel-system solve,
  accepts an optional ``vn`` field on the RHS of the flow-tangency
  rows.
* :func:`aerosim.panel._solve_coupled` — the under-relaxed iteration
  of :eq:`fixed_point`, returning ``Solution.coupled = True`` and
  reporting ``n_iter`` and ``converged``.
* :func:`aerosim.boundary_layer.transpiration_velocity` — equation
  :eq:`vn` with the 1-2-1 smoothing and TE taper applied.
