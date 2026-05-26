The viscous boundary-layer correction
=====================================

Why integral boundary-layer methods?
------------------------------------

The panel method of :doc:`panel_method` solves the *outer*, inviscid
flow. Real airfoils have a thin viscous **boundary layer** hugging the
surface, in which the no-slip condition reduces the velocity from the
inviscid edge value :math:`U_e` to zero at the wall. This layer is the
only place viscosity does anything significant, and it produces

* **skin-friction drag** along the wetted surface;
* **form drag** because the layer thickens the effective body and
  slightly displaces the outer streamlines;
* **separation** and (eventually) **stall** when the adverse pressure
  gradient over the rear of the airfoil overcomes the boundary-layer
  momentum.

A full Navier–Stokes simulation would resolve all this — at the cost of
a 2D mesh, turbulence model, and time-stepping scheme that puts it well
beyond an interactive web demo. **Integral boundary-layer methods** are
a much cheaper compromise: under Prandtl's boundary-layer assumptions
(:math:`\partial p/\partial n \approx 0`, :math:`\delta \ll c`), the
boundary-layer PDE collapses to *ordinary* differential equations along
the surface coordinate :math:`s` for a handful of integral thicknesses.

The integral quantities `aerosim` tracks are

.. math::

   \delta^*(s) &= \int_0^\infty \!\!\Bigl(1 - \tfrac{u}{U_e}\Bigr)\, dn
            \qquad \text{(displacement thickness)}, \\
   \theta(s)   &= \int_0^\infty \!\!\tfrac{u}{U_e}\Bigl(1 - \tfrac{u}{U_e}\Bigr)\, dn
            \qquad \text{(momentum thickness)}, \\
   H(s)        &= \delta^*(s)\,/\,\theta(s)
            \qquad \text{(shape factor)}.

:math:`\delta^*` is how much the outer streamlines are pushed outward by
the layer; it is what eventually couples back into the inviscid solve in
:doc:`coupling`.  :math:`\theta` controls drag — at the trailing edge it
maps directly to the wake momentum deficit and hence to profile drag.
The shape factor :math:`H` is the dimensionless gauge of how close the
layer is to separation: low :math:`H` (~1.4) is healthy turbulent flow,
high :math:`H` (~3–4) is on the edge of separating.

Non-dimensionalisation
~~~~~~~~~~~~~~~~~~~~~~

Everything in the solver is expressed in **chord units** (:math:`c=1`)
and **free-stream speed** units (:math:`V_\infty=1`), so the kinematic
viscosity is

.. math:: \nu = 1/\mathrm{Re},

with the chord-based Reynolds number :math:`\mathrm{Re} = V_\infty c /
\nu`. The edge velocity :math:`U_e(s) = |V_t(s)|/V_\infty` comes
directly from the inviscid solution.

Stagnation point and surface splitting
--------------------------------------

The boundary layer marches **away from the stagnation point** on each
side of the airfoil — upper surface forward to TE, lower surface forward
to TE.  The stagnation point itself is found from the inviscid surface
tangential velocity: it is the point at which :math:`V_t` changes sign,
nearest the leading edge.  Because the node ordering is
TE → upper → LE → lower → TE, this sign change is always near the LE
panels and is unambiguous.

The arc length :math:`s` is then measured from the stagnation point
along each surface separately.


Laminar run — Thwaites' method
------------------------------

For laminar flow with an arbitrary pressure gradient, Thwaites [Thw1949]_
showed empirically that the momentum thickness obeys a simple closed-form
ODE that integrates to

.. math::
   :label: thwaites

   \theta^2(s)\,U_e^6(s)
   \;=\; \theta_0^2 U_{e,0}^6
        \;+\; 0.45\,\nu \int_0^s U_e^5(s')\, ds'.

Starting from :math:`\theta(0) = 0` at the stagnation point — actually,
from the leading-edge similarity value — this is just a cumulative
trapezoidal integral. Given :math:`\theta(s)`, the **local pressure
gradient parameter** is

.. math::
   :label: lambda

   \lambda(s) = \frac{\theta^2}{\nu}\,\frac{dU_e}{ds}.

Thwaites' tabulated correlations then give the local shape factor and
skin friction:

.. math::

   H(\lambda) \approx 2.0\;\dots\;3.55, \qquad
   \ell(\lambda) = \tfrac{\theta\,\tau_w}{\mu\, U_e}.

Thwaites' criterion for laminar **separation** is

.. math::
   :label: thwaites_sep

   \lambda \;\le\; -0.09.

When :eq:`thwaites_sep` is hit the laminar branch is truncated and
turbulent flow is forced from that point onward.


Transition — Michel's criterion
-------------------------------

Real boundary layers do not stay laminar indefinitely; they transition
to turbulence once disturbances grow large enough. `aerosim` uses
**Michel's criterion** [Mic1951]_, the simplest engineering correlation:
transition occurs where

.. math::
   :label: michel

   \mathrm{Re}_\theta \;\ge\;
        1.174\,\bigl(1 + 22\,400/\mathrm{Re}_x\bigr)\,\mathrm{Re}_x^{0.46},

with :math:`\mathrm{Re}_\theta = U_e \theta / \nu` and
:math:`\mathrm{Re}_x = U_e x / \nu`. As :math:`\mathrm{Re}_x` rises with
chord position, the right-hand side approaches :math:`1.174\,
\mathrm{Re}_x^{0.46}`; once the boundary-layer momentum-thickness
Reynolds number catches it, the layer transitions.

The transition location :math:`x_\text{tr}/c` moves **forward as
Reynolds number rises** — a trend the validation suite explicitly
checks.

If laminar separation happens before Michel's criterion fires (as on
heavily loaded leading edges at low :math:`\mathrm{Re}`), the
transition is *forced* at the laminar separation point — a simplified
model for the "laminar separation bubble" that re-attaches downstream
as a turbulent layer.


Turbulent run — Head's entrainment method
-----------------------------------------

Past transition `aerosim` switches to **Head's entrainment method**
[Hea1958]_, which adds a second ODE for an "entrainment shape factor"
:math:`H_1 = (\delta - \delta^*)/\theta` and closes the system with an
empirical entrainment rate :math:`F(H_1)` and a correlation
:math:`H(H_1)`. The two coupled ODEs marched along :math:`s` are

.. math::
   :label: head

   \begin{aligned}
   \frac{d\theta}{ds}     &= \frac{C_f}{2} - (H+2)\,\frac{\theta}{U_e}\,\frac{dU_e}{ds}, \\
   \frac{d(U_e \theta H_1)}{ds} &= U_e\, F(H_1).
   \end{aligned}

The skin-friction closure used here is **Ludwieg–Tillmann** [LT1950]_:

.. math::
   :label: lt

   C_f \;=\; 0.246 \cdot 10^{-0.678\,H}\,\mathrm{Re}_\theta^{-0.268}.

Turbulent **separation** is declared when the shape factor crosses

.. math::

   H \;\gtrsim\; 2.6 \quad (\,\text{equivalently}\;H_1 \;\lesssim\; 3.5\,).

A subtlety, called out in the source code: a closed (cusped) trailing
edge makes :math:`U_e` drop steeply over the last fraction of a percent
of chord, which spikes :math:`H` artificially.  `aerosim` therefore only
flags **genuine separation** when it occurs forward of
:math:`x/c = 0.90`. Trailing-edge separation inside the last 10 % of
chord is treated as a TE-cusp artifact.


Profile drag — the Squire–Young formula
---------------------------------------

The wake far downstream of the trailing edge carries a momentum deficit
that, by the momentum theorem, equals the profile drag. **Squire and
Young** [SY1937]_ derived a closed-form extrapolation that links the
*trailing-edge* boundary-layer state to the *far-wake* momentum
thickness, assuming the wake relaxes to the free stream by a simple
similarity law:

.. math::
   :label: squire_young

   \frac{\theta_\infty}{\theta_\text{TE}}
   \;=\; \left(\frac{U_{e,\text{TE}}}{V_\infty}\right)^{(H_\text{TE} + 5)/2}.

Summed over upper and lower surfaces, the profile-drag coefficient is

.. math::
   :label: cd_visc

   C_{d,\text{visc}}
   \;=\; \frac{2\,\theta_{\infty,\text{up}} + 2\,\theta_{\infty,\text{lo}}}{c}.

The factor of 2 is the chord-based dimensionalisation; with :math:`c=1`
it simply means each surface contributes its far-wake :math:`\theta`.
This is the number stored as ``Solution.cd_visc``, distinct from the
pressure-integrated ``Solution.cd`` (which stays ≈ 0 unless coupling is
on).


Outputs and physical trends
---------------------------

For each surface the solver returns

* :math:`\theta(s)`, :math:`\delta^*(s)`, :math:`H(s)`, :math:`C_f(s)`;
* the **transition location** :math:`x_\text{tr}/c` (marked on the
  :math:`C_p` plot in the web UI);
* a **separation flag** plus the separation arc-length if separated forward of
  the TE cusp;
* the **profile-drag contribution** from Squire–Young.

The validation suite asserts several non-tight but physically essential
trends — these are not coincidences but consequences of the equations
above:

* :math:`C_d` **falls** with increasing :math:`\mathrm{Re}` (thinner
  boundary layer);
* :math:`C_d` **rises** with angle of attack (stronger adverse pressure
  gradient on the upper surface);
* the **transition point moves forward** as :math:`\mathrm{Re}` rises
  (Michel, equation :eq:`michel`);
* skin friction is always **less** than total profile drag
  (form-drag contribution is non-negative).


What this correction does *not* do (yet)
----------------------------------------

The method above is **one-way (uncoupled)**: the inviscid edge velocity
:math:`U_e(s)` is what the panel solver gave us, and the boundary-layer
output never feeds back into it. Lift, surface pressure, and the
d'Alembert :math:`C_d` are therefore *exactly* what the inviscid solve
produced.

That is fine for a drag estimate but misses **viscous decambering** —
the small lift loss caused by the boundary layer effectively thickening
and de-rotating the body. The next chapter, :doc:`coupling`, closes the
loop.
