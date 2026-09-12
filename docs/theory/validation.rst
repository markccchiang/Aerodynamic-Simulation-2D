Validation
==========

The script ``src/validate.py`` is the project's test suite.  It is a
flat list of ``check()`` assertions against analytical and reference
aerodynamics and runs in well under a second.  Each assertion below
corresponds to a single ``check`` in the file; they are grouped here by
the physics they exercise.

Run it after any change to ``panel.py``, ``airfoil.py``,
``boundary_layer.py``, or ``airfoil_io.py``:

.. code-block:: bash

   uv run python src/validate.py


Inviscid solver
---------------

Thin-airfoil and Kutta–Joukowski sanity checks.

* **Symmetric airfoil at zero incidence has zero lift.**
  A NACA 0012 at :math:`\alpha = 0` produces :math:`|C_\ell| < 10^{-3}`.

* **Lift-curve slope.** :math:`dC_\ell/d\alpha` matches the
  thickness-corrected thin-airfoil prediction

  .. math::
     \frac{dC_\ell}{d\alpha} \;\approx\; 2\pi\,(1 + 0.77\,t/c)
                                        \quad\text{per radian},

  to within a few percent.  Note the **thickness correction** — the
  plain :math:`2\pi` flat-plate value is only the zero-thickness limit.

* **Cambered NACA 2412.** At :math:`\alpha = 0` the lift coefficient is
  near 0.25 and the quarter-chord moment is nose-down (negative),
  matching textbook values.

* **d'Alembert's paradox.** The pressure-integrated drag coefficient
  ``Solution.cd`` is essentially zero for any geometry, any angle, with
  the uncoupled BL on or off. This is a *correctness* check, not a bug.

* **Peak surface :math:`C_p \approx 1`.** At the stagnation point
  Bernoulli predicts :math:`C_p = 1`; the discrete solve recovers it.

* **Exact case: flow past a circle.** Potential flow past a circle has
  the closed-form surface pressure

  .. math::
     C_p(\theta) \;=\; 1 - 4\sin^2\theta ,

  and the solver reproduces it to :math:`4.8 \times 10^{-14}` — machine
  precision. Every other assertion in this file is a trend or a
  published number; this one is **analysis**, so it pins
  :func:`~aerosim.panel.induced` and the panel assembly directly rather
  than bounding them. It is built in code at full precision: the same
  comparison against the bundled ``circle.dat`` is limited to
  :math:`4.2 \times 10^{-4}` by the file's six-decimal coordinates, not
  by the solver.

* **Lift from pressure ≡ lift from circulation.**
  :math:`C_\ell` integrated from :math:`C_p` matches the
  Kutta–Joukowski value :math:`C_\ell = -2\,\gamma\,p / V_\infty c`
  (with :math:`p` the panel perimeter and the minus sign reflecting the
  CCW-positive convention) to a few :math:`\times 10^{-3}`.


Uncoupled boundary-layer correction
-----------------------------------

* **Magnitude.** Profile :math:`C_d` for NACA 0012 at
  :math:`\mathrm{Re} = 10^6` lies inside a published band around
  :math:`8 \times 10^{-3}`. This is the only **absolute** tolerance in
  the BL block — everything else is a trend check, because integral
  methods are only good to a few percent.

* **Re trend.** Profile :math:`C_d` falls monotonically as
  :math:`\mathrm{Re}` rises across the usable range.

* **:math:`\alpha` trend.** Profile :math:`C_d` rises with angle of
  attack (stronger adverse pressure gradient on the suction side).

* **Transition advances with Re.** :math:`x_\text{tr}/c` decreases as
  :math:`\mathrm{Re}` rises — a direct consequence of Michel's criterion.

* **Friction < total profile drag.** The skin-friction contribution is
  always less than the Squire–Young total (form drag is non-negative).

* **Inviscid invariance.** Calling ``solve(..., re=Re)`` without
  ``couple=True`` must not change :math:`C_\ell` or the d'Alembert
  :math:`C_d` by **any** detectable amount. This is the contract that
  the uncoupled path is a pure post-process; it is what makes the
  d'Alembert teaching point work.

* **Attached at low :math:`\alpha`, separation onset near stall.** The
  ``separated`` flag is False for healthy operating points and flips on
  as the airfoil approaches stall.


Two-way coupling
----------------

* **Coupling reduces :math:`C_\ell`.** At
  :math:`\alpha = 5°, \mathrm{Re} = 10^6` the coupled
  :math:`C_\ell` is below the inviscid value (viscous decambering).

* **Suction peak softens.** The most-negative :math:`C_p` on the upper
  surface is **less negative** when coupling is on.

* **Lift loss grows as Re falls.** Halving :math:`\mathrm{Re}` widens
  the gap between inviscid and coupled :math:`C_\ell`.

* **Symmetric airfoil at :math:`\alpha = 0` stays at zero lift.**
  Coupling must not break up/down symmetry.

* **Attached case converges.** ``Solution.converged`` is True for a
  case chosen to be well inside the direct-coupling regime.

* **:math:`C_d` is now positive.** Once coupled, the
  pressure-integrated drag has a small non-zero **form-drag**
  component induced by the displacement body — exactly what
  :doc:`coupling` predicts. (The profile drag in
  ``Solution.cd_visc`` is the headline number, but ``Solution.cd`` is
  no longer a d'Alembert zero.)

* **``couple=True`` without ``re=`` raises.** The dispatcher refuses
  to run without a Reynolds number.


.dat file loading
-----------------

* **Selig and Lednicer round-trips.** A NACA 0012 expressed as a Selig
  ``.dat`` file and as a Lednicer ``.dat`` file each load, re-panel,
  and solve to the same :math:`C_\ell` as the geometry built directly
  from the airfoil generator.

* **Normalisation is invariant to scale/offset.** Multiplying all
  coordinates by 5 and shifting by 100 does not change the predicted
  :math:`C_\ell`.

* **The bundled circle is still the analytic case.** Read back from
  ``circle.dat``, :math:`C_p` matches :math:`1 - 4\sin^2\theta` to
  :math:`4.2 \times 10^{-4}` — the file's coordinate precision, which
  makes this a guard on the shipped file rather than on the solver.

* **Every bundled sample loads, re-panels, and solves.** A smoke test
  over all files in ``src/aerosim/airfoils/``.


Web payload and solver plumbing
-------------------------------

The layers *above* the physics, which the rest of the suite does not
touch.

* **The** :math:`C_p` **payload splits at the true leading edge.** The
  upper/lower split uses :math:`\arg\min x` over the control points,
  not the midpoint index. The two coincide only for a symmetric
  section; elsewhere the leading edge sits one or two panels away, and
  splitting at the midpoint would plot lower-surface points on the
  upper curve exactly at the suction peak.

* **Out-of-range requests are rejected by the validation layer.**
  Bad :math:`\alpha`, :math:`\mathrm{Re}` or panel counts become a 422
  from the request model rather than an exception from inside the
  solver.

* **The cambered NACA code snaps** :math:`P: 0 \to 1`. Camber needs a
  non-zero position, so ``M > 0`` with ``P = 0`` cannot silently return
  a symmetric section.

* **A cached geometry matches a fresh solve.** The panel matrix is
  factorized once per :class:`~aerosim.panel.Geometry` and reused (see
  :doc:`panel_method`); a reused geometry must give bit-identical
  results at every angle of attack.

* **Far field → free stream, body masked.** The reconstructed velocity
  field relaxes to :math:`(\cos\alpha, \sin\alpha)` far upstream and
  is NaN inside the body, both of which the web front-end relies on.


Headless web smoke test
-----------------------

Independent of ``validate.py``, the web layer can be exercised without
a browser by importing the FastAPI app and calling the response
builder directly:

.. code-block:: bash

   PYTHONPATH=src uv run python -c "from aerosim.webapp import _pack; \
       from aerosim import naca4, Geometry; \
       print(_pack(Geometry(*naca4('2412')), 'NACA 2412', 5.0, 1e6)['coeffs'])"

If the coefficients block prints, the full solve → BL → JSON pipeline
is healthy end to end.


Bluff bodies: the circle and the golf ball
------------------------------------------

Two of the bundled samples are not airfoils. ``circle.dat`` is the
section of a smooth sphere; ``golfball.dat`` is the *same* circle
carrying 30 raised-cosine dimples, :math:`6 \times 10^{-3}\,c` deep and
:math:`0.082\,c` wide — the proportions of a real ball (a 0.010 in
dimple on a 1.68 in sphere). They are a controlled pair, identical in
every respect except the dimples, and they earn their place for two
opposite reasons.

The first is that the circle is the only body here with a **closed-form
solution**, which is what makes the exact check above possible.

The second is that the pair marks, concretely, where this model stops
being physics. At :math:`\alpha = 0`, :math:`\mathrm{Re} = 10^6`, 200
panels:

.. list-table::
   :header-rows: 1
   :widths: 34 18 18 14 16

   * -
     - :math:`C_d` pressure
     - :math:`C_d` profile
     - :math:`C_\ell`
     - separation
   * - Circle (smooth)
     - 0.00000
     - 0.0006
     - 0.0000
     - :math:`0.893\,c`
   * - Golf ball (dimpled)
     - −0.00000
     - 0.0937
     - 0.0000
     - :math:`0.009\,c`

**Both pressure drags are zero.** That is d'Alembert again
(:doc:`governing_equations`), and it is the whole point of the pair:
dimples cannot change a drag that is identically zero for *any* closed
body in potential flow. Measured values are roughly :math:`C_d = 0.5`
for a smooth sphere and :math:`0.25` for a dimpled one; neither is
recoverable here.

**The separation column runs backwards.** Real dimples *delay*
separation — they trip the boundary layer turbulent, which resists the
adverse gradient and holds on to roughly :math:`115°` from the
stagnation point instead of :math:`82°`, shrinking the wake and halving
the drag. That is the entire reason golf balls have dimples. This model
moves separation the *wrong way*: the smooth circle holds on absurdly
late (:math:`0.893\,c`, about :math:`155°`) and the dimpled one lets go
almost immediately (:math:`0.009\,c`, about :math:`11°`), because each
dimple's steep local adverse gradient trips Thwaites' laminar criterion
on contact. Nothing in an attached-flow integral method can represent
the turbulent reattachment *inside* a dimple that does the real work.
See :doc:`boundary_layer` for why the method is built that way.

**Lift on a circle is an artifact.** A circle has no angle of attack —
rotating it returns the same body — so :math:`C_\ell` must vanish at
every :math:`\alpha`. It does at :math:`\alpha = 0` by symmetry, but the
solver returns :math:`C_\ell = +1.08` at :math:`5°` and :math:`+2.16` at
:math:`10°`. The Kutta condition is anchored to the arbitrary node at
:math:`\theta = 0` and manufactures circulation to put a stagnation
point there. It superficially resembles Magnus lift from backspin; it is
not, and it is the one output of these two samples most likely to be
mistaken for a result.

What the dimples *do* change here is the surface pressure, and that part
is honest potential flow: short-wavelength waviness perturbs the
velocity by roughly :math:`2\pi a/\lambda \approx 36\%`, even though the
amplitude is only :math:`0.6\%` of the diameter.

.. plot::
   :alt: Surface relief and surface pressure for the smooth circle and the dimpled golf ball.
   :caption: The bundled pair at :math:`\alpha = 0`. **Left:** surface
       relief, :math:`r - R`, showing the 30 dimples against the flat
       smooth circle. **Right:** surface pressure. The smooth circle
       (blue) lies exactly on the analytic :math:`1 - 4\sin^2\theta`
       (grey, drawn underneath), while the dimples swing :math:`C_p` by
       nearly :math:`\pm 0.8` about it and drive the shoulder minimum
       from :math:`-3.0` to :math:`-3.8`. Both bodies none the less have
       zero pressure drag.

   import numpy as np
   import matplotlib.pyplot as plt
   from aerosim import Geometry, solve
   from aerosim.airfoil_io import parse_dat, repanel, load_sample_text

   fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.2, 2.9))

   t = np.linspace(0, 360, 721)
   ax2.plot(t, 1 - 4 * np.sin(np.radians(t)) ** 2, color="0.6", lw=2.6,
            label=r"$1-4\sin^2\theta$", zorder=1)

   for key, label, colour in (("circle", "smooth circle", "#1f77b4"),
                              ("golfball", "golf ball", "#d1495b")):
       x, y, _ = parse_dat(load_sample_text(key))
       g = Geometry(*repanel(x, y, 400))
       cp = solve(g, 0.0).cp
       th = np.degrees(np.mod(np.arctan2(g.yc, g.xc - 0.5), 2 * np.pi))
       o = np.argsort(th)
       ax1.plot(th[o], (np.hypot(g.xc - 0.5, g.yc)[o] - 0.5) * 1e3,
                color=colour, lw=1.1, label=label)
       ax2.plot(th[o], cp[o], color=colour, lw=1.0, label=label, zorder=2)

   ax1.set(xlabel=r"$\theta$ (deg)", ylabel=r"$(r-R)\times 10^{3}\ /\ c$",
           xlim=(0, 360), title="surface relief")
   ax2.set(xlabel=r"$\theta$ (deg)", ylabel=r"$C_p$", xlim=(0, 360),
           title=r"surface pressure")
   for a in (ax1, ax2):
       a.set_xticks([0, 90, 180, 270, 360])
   ax2.legend(fontsize=7, frameon=False, loc="lower center")
   fig.tight_layout()

