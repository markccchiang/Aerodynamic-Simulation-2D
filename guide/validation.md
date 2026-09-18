# Validation

← [Back to the README](../README.md)

This is the project's test suite — 37 assertions, **run `uv run python src/validate.py` after any change to the physics core** (`src/aerosim/panel.py`, `airfoil.py`, `boundary_layer.py`, or `airfoil_io.py`); it runs in well under a second.

`validate.py` confirms the solver against analytical / reference results:

- symmetric airfoil produces zero lift at zero incidence;
- lift-curve slope ≈ `2π(1 + 0.77·t/c)` per radian;
- cambered NACA 2412 gives `Cl ≈ 0.25` at α = 0 and a nose-down `Cm`;
- pressure drag ≈ 0 (d'Alembert);
- peak surface `Cp ≈ 1` at the stagnation point;
- lift from the pressure integral matches the Kutta–Joukowski circulation;
- **flow past a circle matches the closed-form `Cp = 1 - 4sin²θ` to 5×10⁻¹⁴** —
  the one assertion here that is analysis rather than a trend or a published
  number, so it pins the influence coefficients directly.

For the viscous correction it also checks:

- profile `Cd` of NACA 0012 ≈ published data at Re = 10⁶;
- `Cd` falls with Reynolds number and rises with angle of attack;
- transition advances toward the leading edge as Reynolds number rises;
- the inviscid `Cl`/`Cd` are untouched by the uncoupled estimate (`couple=False`);
- the flow is attached at α = 0 and separates as it approaches stall.

For the two-way coupling (`couple=True`):

- coupling reduces `Cl` (viscous decambering) and softens the suction peak;
- the lift loss grows as Reynolds number falls (thicker boundary layer);
- a symmetric airfoil at α = 0 gains no spurious lift, and an attached case
  converges; `couple=True` without a Reynolds number raises.

And for `.dat` loading:

- Selig and Lednicer files round-trip to the same `Cl` as the direct geometry;
- coordinates are normalized to unit chord (scale/offset invariant);
- the bundled `circle.dat` still matches the analytic circle to 4×10⁻⁴ — its
  six-decimal coordinates, not the solver, set that bound;
- every bundled sample loads, re-panels, and solves.

And for the layers above the physics:

- the `Cp` payload splits at the true leading edge, not the midpoint index;
- out-of-range API requests are rejected by the validation layer, not by an
  exception from inside the solver;
- a cambered NACA code snaps `P: 0 → 1`, so `M > 0` with `P = 0` cannot quietly
  return a symmetric section;
- a cached geometry gives bit-identical results to a freshly built one;
- the reconstructed flow field relaxes to the free stream far upstream and is
  masked inside the body.

The [validation chapter](../docs/theory/validation.rst) explains what each group is
testing and why the tolerances are what they are.
