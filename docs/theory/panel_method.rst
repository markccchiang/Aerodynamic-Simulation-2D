The Hess–Smith panel method
===========================

Idea
----

Laplace's equation :math:`\nabla^2 \phi = 0` is linear, so any
combination of elementary potential-flow solutions — sources, sinks,
doublets, vortices — is itself a valid flow. The **panel method**
exploits this by *distributing* such elementary singularities on the
surface of the body and choosing their strengths so that the resulting
flow satisfies the boundary conditions of :doc:`governing_equations`.

The version implemented here is the **Hess–Smith** scheme [Hess1967]_,
the original and still the simplest such method for 2D airfoils:

* the airfoil contour is discretised into :math:`N` straight panels;
* each panel carries a **constant-strength source** of unknown strength
  :math:`\sigma_j` (one per panel);
* a **single vortex strength** :math:`\gamma`, the same on every panel,
  supplies the circulation needed for lift.

That makes :math:`N+1` unknowns. They are closed by :math:`N`
flow-tangency equations (one per panel) plus the single Kutta condition,
giving a square :math:`(N+1)\times(N+1)` linear system that is solved
once.


Surface discretisation
----------------------

Given the :math:`N+1` boundary nodes :math:`(x_i, y_i)`, the
:math:`j`-th panel runs from node :math:`j` to node :math:`j+1`. Its
length, midpoint (the **control point**), and orientation are

.. math::

   L_j         &= \sqrt{(x_{j+1}-x_j)^2 + (y_{j+1}-y_j)^2}, \\
   (x_{c,j}, y_{c,j}) &= \tfrac{1}{2}\bigl(x_j + x_{j+1},\; y_j + y_{j+1}\bigr), \\
   \theta_j    &= \arctan2(y_{j+1}-y_j,\; x_{j+1}-x_j),

and the panel-aligned tangent and inward-of-page normal are

.. math::

   \hat{\mathbf{t}}_j &= (\cos\theta_j, \sin\theta_j), \\
   \hat{\mathbf{n}}_j &= (-\sin\theta_j, \cos\theta_j).

.. note::

   `aerosim` uses **cosine spacing** of nodes (clustered near both
   leading and trailing edges) because the surface pressure gradient is
   largest there. The node order is **trailing edge → upper surface →
   leading edge → lower surface → back to trailing edge** — this
   convention is what makes "first and last panels" in the Kutta
   condition correspond to the two trailing-edge panels.

.. plot::
   :alt: NACA 0012 discretized into flat panels with cosine spacing.
   :caption: Cosine-spaced panel discretization of a NACA 0012 (40 panels
       drawn for clarity; the default solve uses ~160). One panel near the
       leading edge is highlighted with its control point (●), tangent
       :math:`\hat{\mathbf{t}}` (blue), and outward normal
       :math:`\hat{\mathbf{n}}` (red). Note how nodes cluster at both the
       leading and trailing edges, where the pressure gradient is largest.

   import numpy as np
   import matplotlib.pyplot as plt

   from aerosim import naca4, Geometry

   x, y = naca4("0012", n_panels=40)
   geom = Geometry(x, y)

   fig, ax = plt.subplots(figsize=(7.5, 2.6))
   ax.plot(geom.x, geom.y, "o-", color="0.3", lw=0.7, ms=2.6)

   j = 7  # a panel a bit aft of the LE on the upper surface
   ax.plot([geom.xa[j], geom.xb[j]], [geom.ya[j], geom.yb[j]],
           color="C2", lw=2.4, zorder=4)
   ax.plot(geom.xc[j], geom.yc[j], "o", color="C2", ms=6, zorder=5)
   scale = 0.07
   ax.annotate("", xy=(geom.xc[j] + scale * geom.tx[j],
                       geom.yc[j] + scale * geom.ty[j]),
               xytext=(geom.xc[j], geom.yc[j]),
               arrowprops=dict(arrowstyle="->", color="C0", lw=1.6))
   ax.annotate("", xy=(geom.xc[j] + scale * geom.noutx[j],
                       geom.yc[j] + scale * geom.nouty[j]),
               xytext=(geom.xc[j], geom.yc[j]),
               arrowprops=dict(arrowstyle="->", color="C3", lw=1.6))
   ax.annotate(r"$\hat{\mathbf{t}}$",
               (geom.xc[j] + 1.1 * scale * geom.tx[j],
                geom.yc[j] + 1.1 * scale * geom.ty[j]),
               color="C0", fontsize=11)
   ax.annotate(r"$\hat{\mathbf{n}}$",
               (geom.xc[j] + 1.1 * scale * geom.noutx[j],
                geom.yc[j] + 1.1 * scale * geom.nouty[j]),
               color="C3", fontsize=11)
   ax.annotate("control point",
               (geom.xc[j], geom.yc[j]),
               xytext=(geom.xc[j] - 0.10, geom.yc[j] + 0.10),
               fontsize=8, color="C2",
               arrowprops=dict(arrowstyle="-", color="C2", lw=0.6))

   ax.set_aspect("equal")
   ax.set_xlim(-0.05, 1.05)
   ax.set_ylim(-0.18, 0.22)
   ax.set_xlabel("x / c")
   ax.set_ylabel("y / c")
   ax.set_title("Panel discretization")
   plt.tight_layout()


Elementary singularities
------------------------

Constant-strength source panel
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A source of strength :math:`\sigma` per unit length, distributed along a
panel of length :math:`L` aligned with the local :math:`x`-axis, induces
the velocity field

.. math::
   :label: source

   \begin{aligned}
   u_p^{(s)}(x,y) &= \frac{\sigma}{4\pi}\,
        \ln\!\frac{(x)^2 + y^2}{(x-L)^2 + y^2}, \\
   w_p^{(s)}(x,y) &= \frac{\sigma}{2\pi}\,
        \bigl[\arctan2(y, x-L) - \arctan2(y, x)\bigr],
   \end{aligned}

where :math:`(x,y)` are coordinates in the panel's **local** frame
(panel along local :math:`+x`, from 0 to :math:`L`).  The subscript
:math:`p` denotes "in panel-local coordinates".

Constant-strength vortex panel
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A vortex sheet of strength :math:`\gamma` per unit length on the same
panel induces

.. math::
   :label: vortex

   \begin{aligned}
   u_p^{(v)}(x,y) &= -\frac{\gamma}{2\pi}\,
        \bigl[\arctan2(y, x-L) - \arctan2(y, x)\bigr], \\
   w_p^{(v)}(x,y) &= \phantom{-}\frac{\gamma}{4\pi}\,
        \ln\!\frac{(x)^2 + y^2}{(x-L)^2 + y^2}.
   \end{aligned}

Comparing :eq:`source` and :eq:`vortex`, a vortex panel's induced field
is just a source panel's field rotated by 90°.  In `aerosim` the same
function :func:`aerosim.panel.induced` returns both kinds in one call,
which is why the implementation is so compact.

Returning to the global frame is a rotation by :math:`\theta_j`:

.. math::

   \begin{pmatrix} u \\ v \end{pmatrix}
   = \begin{pmatrix}\cos\theta_j & -\sin\theta_j \\
                    \sin\theta_j & \phantom{-}\cos\theta_j\end{pmatrix}
     \begin{pmatrix} u_p \\ w_p \end{pmatrix}.

The two elementary fields look as follows. A **source panel** (left) blows
flow outward symmetrically; a **vortex panel** (right) drives a circulation
around it. The vortex's pattern is the source's pattern rotated by 90°,
which is why one influence routine returns both at once.

.. plot::
   :alt: Velocity fields induced by a unit-strength source panel and vortex panel.
   :caption: Velocity field around a single unit-strength source panel (left)
       and a single unit-strength vortex panel (right). The panel itself runs
       from :math:`(0,0)` to :math:`(1,0)`. Streamlines are coloured by
       speed; the panel is drawn as a thick black bar.

   import numpy as np
   import matplotlib.pyplot as plt

   from aerosim import Geometry
   from aerosim.panel import induced

   # A degenerate one-panel "geometry" just to exercise induced().
   geom = Geometry(np.array([0.0, 1.0]), np.array([0.0, 0.0]))

   xs = np.linspace(-0.6, 1.6, 60)
   ys = np.linspace(-0.8, 0.8, 50)
   X, Y = np.meshgrid(xs, ys)
   us_x, us_y, uv_x, uv_y = induced(X.ravel(), Y.ravel(), geom)

   Us = us_x[:, 0].reshape(X.shape)
   Vs = us_y[:, 0].reshape(X.shape)
   Uv = uv_x[:, 0].reshape(X.shape)
   Vv = uv_y[:, 0].reshape(X.shape)

   fig, (axs, axv) = plt.subplots(1, 2, figsize=(8.4, 3.3))
   axs.streamplot(X, Y, Us, Vs, density=1.2, color=np.hypot(Us, Vs),
                  cmap="magma", linewidth=0.7, arrowsize=0.7)
   axs.plot([0, 1], [0, 0], "-", color="k", lw=3)
   axs.set_aspect("equal"); axs.set_xlabel("x"); axs.set_ylabel("y")
   axs.set_title("Source panel (σ = 1)")

   axv.streamplot(X, Y, Uv, Vv, density=1.2, color=np.hypot(Uv, Vv),
                  cmap="magma", linewidth=0.7, arrowsize=0.7)
   axv.plot([0, 1], [0, 0], "-", color="k", lw=3)
   axv.set_aspect("equal"); axv.set_xlabel("x"); axv.set_ylabel("y")
   axv.set_title("Vortex panel (γ = 1)")

   plt.tight_layout()


Influence coefficients
----------------------

Let :math:`(u_{ij}^{(s)},\, v_{ij}^{(s)})` denote the velocity induced
at control point :math:`i` by a **unit-strength** source on panel
:math:`j`, and likewise :math:`(u_{ij}^{(v)},\, v_{ij}^{(v)})` for a
unit vortex. These four :math:`N\times N` matrices are the *influence
coefficients* — they encode the geometry once and for all.

On the diagonal (a panel inducing flow at *its own* control point) the
self-influence is obtained as the limit from the **flow side** of the
panel. For a source it is a half-strength outward-normal velocity; for
a vortex, a half-strength tangential velocity:

.. math::

   \mathbf{V}_{ii}^{(s)} = \tfrac{1}{2}\sigma\, \hat{\mathbf{n}}_i^\text{ext},
   \qquad
   \mathbf{V}_{ii}^{(v)} = -\tfrac{1}{2}\gamma\, \hat{\mathbf{t}}_i.

Which side is the "flow side" depends on the polygon orientation — for
the TE→upper→LE→lower→TE ordering used here, it is detected
automatically in :meth:`aerosim.panel.Geometry._detect_exterior_side` by
a point-in-polygon test.


The linear system
-----------------

The total velocity at control point :math:`i` is the free stream
:math:`\mathbf{V}_\infty = V_\infty(\cos\alpha, \sin\alpha)` plus all
panel contributions:

.. math::
   :label: total

   \mathbf{V}_i \;=\; \mathbf{V}_\infty
                  \,+\, \sum_{j=1}^{N} \sigma_j\, \mathbf{V}_{ij}^{(s)}
                  \,+\, \gamma \sum_{j=1}^{N} \mathbf{V}_{ij}^{(v)}.

Flow tangency  (N equations)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For each control point :math:`i`, equation
:math:`\mathbf{V}_i \cdot \hat{\mathbf{n}}_i = 0` becomes

.. math::
   :label: tang_lin

   \sum_{j=1}^{N} \sigma_j\,(\mathbf{V}_{ij}^{(s)} \cdot \hat{\mathbf{n}}_i)
   \,+\, \gamma \sum_{j=1}^{N}(\mathbf{V}_{ij}^{(v)} \cdot \hat{\mathbf{n}}_i)
   \;=\; -\, \mathbf{V}_\infty \cdot \hat{\mathbf{n}}_i.

That is :math:`N` linear equations for :math:`N+1` unknowns.

Kutta condition  (1 equation)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The remaining equation enforces that the upper- and lower-surface
tangential velocities at the trailing edge are equal in magnitude and
opposite in sign. In the discrete setting, the cleanest way to express
this is to require that the **sum of the tangential velocities on the
two TE-adjacent panels is zero**:

.. math::
   :label: kutta_lin

   \mathbf{V}_1 \cdot \hat{\mathbf{t}}_1
   \,+\, \mathbf{V}_N \cdot \hat{\mathbf{t}}_N
   \;=\; 0.

Stacking :eq:`tang_lin` for :math:`i=1,\dots,N` and :eq:`kutta_lin` gives
a dense :math:`(N+1)\times(N+1)` linear system

.. math::

   \mathbf{A}\,
   \begin{pmatrix} \sigma_1 \\ \vdots \\ \sigma_N \\ \gamma \end{pmatrix}
   \;=\; \mathbf{b},

solved once with a direct factorisation (`numpy.linalg.solve`).


Post-processing
---------------

Surface velocity and pressure
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

With :math:`(\sigma_j)` and :math:`\gamma` known, evaluate
:eq:`total` at each control point and take the **tangential**
component:

.. math::

   V_{t,i} = \mathbf{V}_i \cdot \hat{\mathbf{t}}_i, \qquad
   C_{p,i} = 1 - \left(\frac{V_{t,i}}{V_\infty}\right)^2.

The normal component is zero by construction.

Force and moment coefficients
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Integrating the pressure (equation :eq:`forces` of
:doc:`governing_equations`) over the discrete contour reduces to a
weighted sum over panels:

.. math::

   \begin{aligned}
   F_x &= -\sum_j C_{p,j}\, n^\text{ext}_{x,j}\, L_j, \\
   F_y &= -\sum_j C_{p,j}\, n^\text{ext}_{y,j}\, L_j, \\
   C_\ell &= -F_x \sin\alpha + F_y \cos\alpha, \\
   C_d    &= \phantom{-}F_x \cos\alpha + F_y \sin\alpha
                  \;\;\approx\; 0
                  \quad \text{(d'Alembert; numerical check)}, \\
   C_{m,c/4} &= \sum_j C_{p,j}\, L_j
                  \bigl[(x_{c,j}-\tfrac{1}{4})\, n^\text{ext}_{y,j}
                        - y_{c,j}\, n^\text{ext}_{x,j}\bigr].
   \end{aligned}

Lift-curve slope and thickness
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Thin-airfoil theory predicts a lift-curve slope of :math:`2\pi` per
radian. For a finite-thickness airfoil a useful refinement (and the one
the validation suite uses) is

.. math::

   \frac{dC_\ell}{d\alpha}
       \;\approx\; 2\pi\,(1 + 0.77\,t/c),

where :math:`t/c` is the maximum thickness-to-chord ratio. A NACA 0012
section therefore comes out near :math:`6.77\;\text{rad}^{-1}`, slightly
above the flat-plate value.


Off-body velocity field
~~~~~~~~~~~~~~~~~~~~~~~

The same equation :eq:`total` evaluated at arbitrary field points
:math:`(x_p, y_p)` gives the velocity anywhere outside the body. That
is what :mod:`aerosim.flowfield` does to build the streamlines drawn in
the web UI — there is no separate "streamline solver"; the panel
strengths plus the free stream *are* the field.


From inviscid to real-airfoil performance
-----------------------------------------

The pipeline so far gives an excellent prediction of :math:`C_\ell` and
:math:`C_p` at small angles of attack, but no drag and no stall: the
limitations of pure potential flow listed at the end of
:doc:`governing_equations`. The next chapter, :doc:`boundary_layer`,
develops the viscous correction that supplies the missing profile drag,
locates transition, and flags incipient separation.
