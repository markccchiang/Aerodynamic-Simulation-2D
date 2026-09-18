# References

← [Back to the README](../README.md)

Where each method in `aerosim` comes from. DOIs link to the publisher's record.
The Sphinx theory docs have a References chapter too, next to the derivations
these sources support (see [Building the Sphinx docs](../README.md#building-the-sphinx-docs)).

## Methods, where they live, and their sources

| Method | Used in | Source |
|--------|---------|--------|
| Hess–Smith source + vortex panel method | [`panel.py`](../src/aerosim/panel.py) | Hess & Smith (1967) |
| Kutta condition; Kutta–Joukowski lift check | [`panel.py`](../src/aerosim/panel.py), [`validate.py`](../src/validate.py) | Kutta (1902); Joukowski (1906) |
| NACA 4-digit thickness and camber lines | [`airfoil.py`](../src/aerosim/airfoil.py) | Abbott & von Doenhoff (1959) |
| Selig and Lednicer `.dat` layouts | [`airfoil_io.py`](../src/aerosim/airfoil_io.py) | UIUC Airfoil Coordinates Database |
| Thwaites' laminar boundary layer | [`boundary_layer.py`](../src/aerosim/boundary_layer.py) | Thwaites (1949) |
| Michel's transition criterion | [`boundary_layer.py`](../src/aerosim/boundary_layer.py) | Michel (1951) |
| Head's entrainment method (turbulent) | [`boundary_layer.py`](../src/aerosim/boundary_layer.py) | Head (1958) |
| Ludwieg–Tillmann skin friction | [`boundary_layer.py`](../src/aerosim/boundary_layer.py) | Ludwieg & Tillmann (1950) |
| Squire–Young profile drag | [`boundary_layer.py`](../src/aerosim/boundary_layer.py) | Squire & Young (1937) |

## Textbooks

Any one of these is enough to follow the method; together they cover everything
from Laplace's equation to viscous–inviscid coupling.

- **Anderson, J. D.** *Fundamentals of Aerodynamics*, 6th ed. McGraw-Hill, 2017.
  Potential flow, thin-airfoil theory and boundary layers at textbook depth.
- **Katz, J., and Plotkin, A.** *Low-Speed Aerodynamics*, 2nd ed. Cambridge
  University Press, 2001. The standard reference for panel methods, including
  constant-strength source and vortex panels.
- **Schlichting, H., and Gersten, K.** *Boundary-Layer Theory*, 9th ed. Springer,
  2017. Differential and integral boundary-layer methods, including Thwaites'
  and Head's.
- **Drela, M.** *Flight Vehicle Aerodynamics*. MIT Press, 2014. A modern treatment
  by XFOIL's author, including viscous–inviscid interaction.

## Original method papers

- **Hess, J. L., and Smith, A. M. O.** "Calculation of Potential Flow about
  Arbitrary Bodies." *Progress in Aerospace Sciences*, Vol. 8, 1967, pp. 1–138.
  [doi:10.1016/0376-0421(67)90003-6](https://doi.org/10.1016/0376-0421%2867%2990003-6)
- **Thwaites, B.** "Approximate Calculation of the Laminar Boundary Layer."
  *Aeronautical Quarterly*, Vol. 1, 1949, pp. 245–280.
- **Michel, R.** "Étude de la transition sur les profils d'aile; établissement
  d'un critère de détermination du point de transition et calcul de la traînée
  de profil incompressible." ONERA Rapport 1/1578-A, 1951.
- **Head, M. R.** "Entrainment in the Turbulent Boundary Layer." Aeronautical
  Research Council R&M 3152, 1958.
- **Ludwieg, H., and Tillmann, W.** "Investigations of the Wall-Shearing Stress
  in Turbulent Boundary Layers." NACA TM 1285, 1950 — English translation of the
  1949 German original.
- **Squire, H. B., and Young, A. D.** "The Calculation of the Profile Drag of
  Aerofoils." Aeronautical Research Council R&M 1838, 1937.
- **Goldstein, S.** "On Laminar Boundary-Layer Flow Near a Position of
  Separation." *Quarterly Journal of Mechanics and Applied Mathematics*, Vol. 1,
  1948, pp. 43–69. [doi:10.1093/qjmam/1.1.43](https://doi.org/10.1093/qjmam/1.1.43)
  The singularity at separation that stops a direct viscous–inviscid coupling
  from marching through it — the limit described in [How it works](how-it-works.md).

## Airfoil data and file formats

- **Abbott, I. H., and von Doenhoff, A. E.** *Theory of Wing Sections: Including
  a Summary of Airfoil Data*. Dover, 1959. The NACA 4-digit thickness and camber
  equations that `naca4()` implements, and tabulated section data for checking
  predictions against.
- **UIUC Airfoil Coordinates Database** (Selig, University of Illinois).
  [m-selig.ae.illinois.edu/ads/coord_database.html](https://m-selig.ae.illinois.edu/ads/coord_database.html)
  Over a thousand airfoils as `.dat` files, most in the Selig layout that
  [Loading airfoils](airfoils.md) describes.

## Bluff bodies

Sources for the real-world behaviour that the
[golf-ball section](airfoils.md#a-bluff-body-the-bundled-golf-ball) contrasts
with what the model can show.

- **Achenbach, E.** "Experiments on the Flow Past Spheres at Very High Reynolds
  Numbers." *Journal of Fluid Mechanics*, Vol. 54, 1972, pp. 565–575.
  [doi:10.1017/S0022112072000874](https://doi.org/10.1017/S0022112072000874)
  Sphere drag and separation angle through the drag crisis.
- **Bearman, P. W., and Harvey, J. K.** "Golf Ball Aerodynamics." *Aeronautical
  Quarterly*, Vol. 27, 1976, pp. 112–122.
  [doi:10.1017/S0001925900007617](https://doi.org/10.1017/S0001925900007617)
  Drag and lift measurements on dimpled balls.

## Historical background

- **d'Alembert, J. le R.** *Essai d'une nouvelle théorie de la résistance des
  fluides*, 1752. The original statement of the paradox that the panel method
  reproduces and the boundary-layer correction repairs.
- **Kutta, M. W.** "Auftriebskräfte in strömenden Flüssigkeiten." *Illustrierte
  Aeronautische Mitteilungen*, 1902.
- **Joukowski, N. E.** "Sur les tourbillons adjoints." *Trans. Phys. Sect. Imp.
  Soc. Friends of Nat. Sci. (Moscow)*, 1906. With Kutta's paper, the
  circulation–lift theorem `L' = ρ V∞ Γ`.
- **Prandtl, L.** "Über Flüssigkeitsbewegung bei sehr kleiner Reibung."
  *Verhandlungen des III. Internationalen Mathematiker-Kongresses*, Heidelberg,
  1904. The original boundary-layer paper.

## Software for comparison

- **XFOIL** (Drela). [web.mit.edu/drela/Public/web/xfoil](https://web.mit.edu/drela/Public/web/xfoil/)
  The reference 2D airfoil code: a linear-vorticity panel method coupled to an
  integral boundary layer with e<sup>N</sup> transition, with the viscous and
  inviscid equations solved simultaneously by Newton's method. `aerosim`
  implements a deliberately simpler subset of the same ideas. Described in
  **Drela, M.** "XFOIL: An Analysis and Design System for Low Reynolds Number
  Airfoils." In *Low Reynolds Number Aerodynamics*, Lecture Notes in Engineering,
  Vol. 54, Springer, 1989, pp. 1–12.
  [doi:10.1007/978-3-642-84010-4_1](https://doi.org/10.1007/978-3-642-84010-4_1)
