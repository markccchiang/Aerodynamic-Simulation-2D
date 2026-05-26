Symbol glossary
===============

This is the notation used throughout the theory chapters.  Where the
code uses a different name it is noted in parentheses.

Geometry and free stream
------------------------

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Symbol
     - Meaning
   * - :math:`c`
     - Chord length. Non-dimensionalisation reference; the code uses
       :math:`c = 1` throughout.
   * - :math:`t/c`
     - Maximum thickness-to-chord ratio.
   * - :math:`V_\infty`
     - Free-stream speed (``vinf`` in code; default 1).
   * - :math:`\alpha`
     - Geometric angle of attack (``alpha_deg`` in code, degrees).
   * - :math:`\mathrm{Re}`
     - Chord-based Reynolds number, :math:`V_\infty c / \nu`.
   * - :math:`\nu`
     - Kinematic viscosity. In code units, :math:`\nu = 1/\mathrm{Re}`.
   * - :math:`(x, y)`
     - Cartesian airfoil coordinates, unit chord.
   * - :math:`s`
     - Surface arc length (measured from the stagnation point inside
       the boundary-layer routines).

Panel discretisation
--------------------

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Symbol
     - Meaning
   * - :math:`N`
     - Number of panels.
   * - :math:`(x_{c,j}, y_{c,j})`
     - Control point of panel :math:`j` (panel midpoint).
   * - :math:`L_j`
     - Length of panel :math:`j`.
   * - :math:`\theta_j`
     - Panel orientation angle in the global frame.
   * - :math:`\hat{\mathbf{t}}_j, \hat{\mathbf{n}}_j`
     - Panel-aligned unit tangent and normal.
   * - :math:`\sigma_j`
     - Source strength per unit length on panel :math:`j`
       (``sigma`` in code).
   * - :math:`\gamma`
     - Single vortex strength shared by all panels
       (``gamma`` in code).

Inviscid output fields
----------------------

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Symbol
     - Meaning
   * - :math:`\phi`
     - Velocity potential, :math:`\mathbf{V} = \nabla \phi`.
   * - :math:`V_t(s)`
     - Surface tangential velocity (``vt`` in code).
   * - :math:`C_p(s)`
     - Surface pressure coefficient.
   * - :math:`C_\ell`
     - Lift coefficient (``cl``).
   * - :math:`C_d`
     - Pressure-integrated drag coefficient (``cd``); ≈ 0 unless
       coupling is on (d'Alembert's paradox).
   * - :math:`C_{m,c/4}`
     - Quarter-chord pitching-moment coefficient (``cm_qc``),
       positive = nose-up.
   * - :math:`\Gamma`
     - Total circulation around the airfoil; lift per unit span is
       :math:`\rho V_\infty \Gamma` (Kutta–Joukowski).

Boundary layer
--------------

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Symbol
     - Meaning
   * - :math:`U_e(s)`
     - Inviscid edge velocity at the boundary-layer edge
       (``ue`` in code).
   * - :math:`\theta(s)`
     - Momentum thickness.
   * - :math:`\delta^*(s)`
     - Displacement thickness.
   * - :math:`H(s)`
     - Shape factor, :math:`\delta^*/\theta`.
   * - :math:`H_1(s)`
     - Head's entrainment shape factor, :math:`(\delta - \delta^*)/\theta`.
   * - :math:`C_f(s)`
     - Local skin-friction coefficient, referenced to :math:`U_e`.
   * - :math:`\mathrm{Re}_\theta`
     - Momentum-thickness Reynolds number, :math:`U_e \theta / \nu`.
   * - :math:`\mathrm{Re}_x`
     - Local chord-Reynolds number, :math:`U_e x / \nu`.
   * - :math:`\lambda`
     - Thwaites' pressure-gradient parameter,
       :math:`(\theta^2/\nu)\,dU_e/ds`.
   * - :math:`x_\text{tr}/c`
     - Transition location.
   * - :math:`C_{d,\text{visc}}`
     - Profile drag from Squire–Young (``cd_visc`` / ``bl.cd``).

Coupling
--------

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Symbol
     - Meaning
   * - :math:`v_n(s)`
     - Wall-transpiration velocity, :math:`d(U_e \delta^*)/ds`.
   * - :math:`\omega`
     - Under-relaxation factor in the direct iteration
       (``_COUPLE_RELAX = 0.25``).
   * - :math:`\varepsilon`
     - Coupling tolerance on :math:`\Delta C_\ell`
       (``_COUPLE_TOL = 2\times10^{-5}``).
   * - :math:`k_\text{max}`
     - Coupling iteration cap (``_COUPLE_MAX_ITER = 80``).
   * - :math:`x_\text{taper}`
     - Trailing-edge taper start for :math:`v_n` (``_TE_TAPER_X = 0.95``).
