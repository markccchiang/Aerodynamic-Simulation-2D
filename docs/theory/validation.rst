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

* **Every bundled sample loads, re-panels, and solves.** A smoke test
  over all files in ``src/aerosim/airfoils/``.


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
