References
==========

Core textbooks
--------------

These are the standard references for the theory the simulator
implements.  Any one of them is enough to follow the chapters in this
document; together they cover everything from Laplace's equation to the
viscous–inviscid coupling.

* **Anderson, J. D.** *Fundamentals of Aerodynamics* (6th ed.,
  McGraw-Hill, 2017).
  General textbook treatment of potential flow, thin-airfoil theory,
  boundary layers, and a chapter on integral methods.

* **Katz, J., and Plotkin, A.** *Low-Speed Aerodynamics* (2nd ed.,
  Cambridge University Press, 2001).
  The canonical reference for panel methods of every flavour, with
  detailed derivations of the source / vortex / doublet building
  blocks and the Hess–Smith scheme in particular.

* **Schlichting, H., and Gersten, K.**
  *Boundary-Layer Theory* (9th ed., Springer, 2017).
  The standard reference for both differential and integral
  boundary-layer methods, including Thwaites, Head, and the
  Ludwieg–Tillmann correlation.

* **Drela, M.** *Flight Vehicle Aerodynamics* (MIT Press, 2014).
  Modern treatment that includes wall-transpiration coupling and the
  semi-inverse schemes used by XFOIL.

Original method papers
----------------------

The integral boundary-layer chain implemented in
:mod:`aerosim.boundary_layer` is the classic "XFOIL-lite" sequence, with
each link tracing to a single foundational paper:

.. [Hess1967] Hess, J. L., and Smith, A. M. O.
   "Calculation of Potential Flow about Arbitrary Bodies."
   *Progress in Aerospace Sciences*, Vol. 8, 1967, pp. 1–138.

.. [Thw1949] Thwaites, B.
   "Approximate Calculation of the Laminar Boundary Layer."
   *Aeronautical Quarterly*, Vol. 1, 1949, pp. 245–280.

.. [Mic1951] Michel, R.
   "Étude de la transition sur les profils d'aile;
   établissement d'un critère de détermination du point de transition
   et calcul de la traînée de profil incompressible."
   ONERA Rapport 1/1578-A, 1951.

.. [Hea1958] Head, M. R.
   "Entrainment in the Turbulent Boundary Layer."
   Aeronautical Research Council R&M 3152, 1958.

.. [LT1950] Ludwieg, H., and Tillmann, W.
   "Investigations of the Wall-Shearing Stress in Turbulent Boundary
   Layers."  NACA TM 1285, 1950 (English translation of the 1949
   German original).

.. [SY1937] Squire, H. B., and Young, A. D.
   "The Calculation of the Profile Drag of Aerofoils."
   Aeronautical Research Council R&M 1838, 1937.

.. [Gol1948] Goldstein, S.
   "On Laminar Boundary-Layer Flow Near a Position of Separation."
   *Quarterly Journal of Mechanics and Applied Mathematics*, Vol. 1,
   No. 1, 1948, pp. 43–69.

Background
----------

* **d'Alembert, J. le R.**
  *Essai d'une nouvelle théorie de la résistance des fluides*, 1752.
  The original statement of the paradox that the present panel method
  reproduces and the boundary-layer correction repairs.

* **Prandtl, L.**
  "Über Flüssigkeitsbewegung bei sehr kleiner Reibung."
  *Verhandlungen des III. Internationalen Mathematiker-Kongresses*,
  Heidelberg, 1904. The original boundary-layer paper.

* **Kutta, M. W.**
  "Auftriebskräfte in strömenden Flüssigkeiten."
  *Illustrierte Aeronautische Mitteilungen*, 1902.

* **Joukowski, N. E.**
  "Sur les tourbillons adjoints."
  *Trans. Phys. Sect. Imp. Soc. Friends of Nat. Sci. (Moscow)*,
  1906. Together with Kutta's 1902 paper this gives the
  circulation–lift theorem :math:`L' = \rho V_\infty \Gamma`.

Software for comparison
-----------------------

* **XFOIL** (Drela). The reference 2D airfoil code that combines a
  higher-order panel method with a full viscous–inviscid coupling via
  an :math:`e^N` transition model and Veldman's semi-inverse scheme.
  `aerosim` implements a deliberately simpler subset of the same ideas.

* **NACA airfoil sections** — published tabulated data for
  4-digit, 5-digit, and 6-series airfoils. Useful for sanity-checking
  drag predictions; the figure-of-merit in ``validate.py`` for NACA
  0012 :math:`C_d` at :math:`\mathrm{Re} = 10^6` comes from this
  database.
