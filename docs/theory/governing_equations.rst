Governing equations
===================

Modelling assumptions
---------------------

The flow around the airfoil is taken to be

* **two-dimensional** — variations along the span are ignored;
* **steady** — :math:`\partial/\partial t \equiv 0`;
* **incompressible** — :math:`\rho` is constant, valid for free-stream
  Mach numbers below ~0.3;
* **inviscid** in the bulk — viscosity is confined to a thin boundary
  layer (treated separately in :doc:`boundary_layer`);
* **irrotational** in the bulk — the upstream flow has no vorticity, and
  for an inviscid fluid Kelvin's theorem then keeps the flow irrotational
  everywhere outside thin shear layers.

These assumptions reduce the Navier–Stokes equations to a single linear
PDE: a setting in which a panel method is both possible and natural.


From the conservation laws
--------------------------

Mass conservation for an incompressible fluid is the **continuity
equation**

.. math::
   :label: continuity

   \nabla \cdot \mathbf{V} = 0,

with :math:`\mathbf{V}=(u,v)` the velocity field.  Irrotationality is

.. math::
   :label: irrot

   \nabla \times \mathbf{V} = \mathbf{0}.

Equation :eq:`irrot` is the integrability condition that lets us write
the velocity as the gradient of a scalar **velocity potential**
:math:`\phi`:

.. math::
   :label: potential

   \mathbf{V} = \nabla \phi.

Substituting :eq:`potential` into the continuity equation
:eq:`continuity` yields **Laplace's equation**:

.. math::
   :label: laplace

   \nabla^2 \phi = 0.

That is the fundamental equation `aerosim` solves.  It is linear and
elliptic; superposition of any two solutions is again a solution. This
is what makes the panel method possible: each elementary singularity is
a closed-form solution of :eq:`laplace`, and the body's contribution can
be written as a *sum* of such elementary fields.


Boundary conditions
-------------------

To pick out the unique flow around our airfoil we attach three boundary
conditions:

1. **Free-stream at infinity.** Far from the body the flow tends to the
   uniform velocity

   .. math::
      \mathbf{V} \;\to\; \mathbf{V}_\infty
                    = (V_\infty \cos\alpha,\;V_\infty \sin\alpha)
      \qquad\text{as}\quad |\mathbf{x}| \to \infty,

   where |alpha| is the geometric angle of attack.

2. **Flow tangency on the body.** A solid surface cannot have flow
   passing through it. With :math:`\hat{\mathbf{n}}` the unit outward
   normal,

   .. math::
      :label: tangency

      \mathbf{V} \cdot \hat{\mathbf{n}} = 0 \quad \text{on the airfoil}.

   Equation :eq:`tangency` discretises into one equation per panel in
   :doc:`panel_method`.

3. **The Kutta condition.** Inviscid flow around a sharp-trailing-edge
   airfoil admits a one-parameter family of solutions, parameterised by
   the circulation :math:`\Gamma = \oint \mathbf{V}\cdot d\mathbf{l}`.
   The physical choice — the one that matches real, viscous flow at
   small angles — is the one for which **the flow leaves the trailing
   edge smoothly**: no infinite velocity at the sharp corner, and the
   upper- and lower-surface tangential velocities at the trailing edge
   are equal in magnitude and opposite in sign.

   This **Kutta condition** picks out a unique :math:`\Gamma`, and by
   the **Kutta–Joukowski theorem** the lift per unit span is

   .. math::
      :label: kj

      L' = \rho_\infty V_\infty \Gamma.


The pressure coefficient
------------------------

Once :math:`\phi` (and hence :math:`\mathbf{V}`) is known, the surface pressure
follows from **Bernoulli's equation**, which for a steady, inviscid,
incompressible flow reads

.. math::
   p + \tfrac{1}{2}\rho |\mathbf{V}|^2 = p_\infty + \tfrac{1}{2}\rho V_\infty^2.

The dimensionless **pressure coefficient** is then

.. math::
   :label: cp

   C_p \;=\; \frac{p - p_\infty}{\tfrac{1}{2}\rho V_\infty^2}
      \;=\; 1 - \left(\frac{|\mathbf{V}|}{V_\infty}\right)^2.

On the body, where the only surviving velocity component is the surface
tangential velocity :math:`V_t`, this collapses to
:math:`C_p = 1 - (V_t / V_\infty)^2`, which is exactly what
:func:`aerosim.panel._solve_panels` evaluates at every control point.

A representative inviscid solution from `aerosim` is shown below: the
streamline pattern around a NACA 2412 at :math:`\alpha = 5°` and the
corresponding surface :math:`C_p` distribution. The upper surface
(lower curve) shows a strong suction peak near the leading edge; the
lower surface (upper curve) carries higher pressure. Their integrated
difference is the lift.

.. plot::
   :alt: Streamlines and surface pressure coefficient for NACA 2412 at alpha=5 deg.
   :caption: Inviscid potential flow around a NACA 2412 at α = 5°: streamlines
       on the left, surface :math:`C_p` (negative upward, aerodynamic
       convention) on the right.

   import numpy as np
   import matplotlib.pyplot as plt

   from aerosim import naca4, Geometry, solve, velocity_field

   geom = Geometry(*naca4("2412"))
   sol = solve(geom, alpha_deg=5.0)

   xs = np.linspace(-0.5, 1.5, 280)
   ys = np.linspace(-0.6, 0.6, 180)
   X, Y, U, V = velocity_field(sol, xs, ys)

   fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8.4, 3.4))

   speed = np.hypot(U, V)
   ax1.streamplot(X, Y, U, V, density=1.6, color=speed, cmap="viridis",
                  linewidth=0.7, arrowsize=0.7)
   ax1.fill(geom.x, geom.y, color="0.15", zorder=5)
   ax1.set_xlim(-0.4, 1.4)
   ax1.set_ylim(-0.5, 0.5)
   ax1.set_aspect("equal")
   ax1.set_xlabel("x / c")
   ax1.set_ylabel("y / c")
   ax1.set_title("Streamlines (NACA 2412, α = 5°)")

   ax2.plot(geom.xc, sol.cp, lw=1.2)
   ax2.axhline(0, color="0.5", lw=0.6)
   ax2.invert_yaxis()
   ax2.set_xlabel("x / c")
   ax2.set_ylabel(r"$C_p$")
   ax2.set_title(f"Surface $C_p$  ($C_\\ell$ = {sol.cl:.3f})")

   plt.tight_layout()

The aerodynamic coefficients follow by integrating the pressure
distribution around the contour:

.. math::
   :label: forces

   \begin{aligned}
   \mathbf{F}\;/\; (\tfrac{1}{2}\rho V_\infty^2 c)
     &= -\oint C_p \,\hat{\mathbf{n}}\, ds, \\
   C_\ell &= -F_x \sin\alpha + F_y \cos\alpha, \\
   C_d    &= \phantom{-}F_x \cos\alpha + F_y \sin\alpha, \\
   C_m    &= \oint C_p \big[(x - x_\text{ref})\, n_y
                            - y\, n_x\big]\, ds.
   \end{aligned}

Here :math:`x_\text{ref} = c/4` is used for the **quarter-chord moment**
:math:`C_{m,c/4}`, and the sign convention is the aerodynamic one
(positive moment = nose-up). With a unit chord :math:`c=1`, this maps
directly to the surface integral computed in the solver.


d'Alembert's paradox
--------------------

A celebrated consequence of solving :eq:`laplace` subject to
:eq:`tangency` and a uniform free stream is that **a closed body in a
steady, inviscid, incompressible, irrotational flow experiences no
drag**. The pressure forces on the windward and leeward sides cancel
exactly. The lift survives because the Kutta condition picks a
*circulatory* component that breaks the up–down symmetry, but no
streamwise force ever appears.

This result was a major puzzle when Jean le Rond d'Alembert published it
in 1752: real airfoils plainly do experience drag.  The resolution is
that *no flow is ever truly inviscid* — there is always a thin viscous
boundary layer near the surface, and downstream of separation a wake.
Boundary-layer theory (Prandtl, 1904) is the bridge between the inviscid
outer flow and reality.

For `aerosim` the practical consequence is:

* The pressure-integrated :math:`C_d` returned by the panel solver is
  **expected to be ~0**. The CLAUDE.md warning that a near-zero
  ``Solution.cd`` is "a correctness check, not a bug" is just
  d'Alembert's paradox in code form.
* The realistic profile-drag estimate must come from elsewhere — the
  integral boundary-layer method developed in :doc:`boundary_layer`.

The next chapter, :doc:`panel_method`, shows how Laplace's equation is
turned into a finite linear system using surface singularities, with the
Kutta condition closing the system to give a single physical answer.
