"""Validate the panel solver against known analytical / reference results.

Run with:  uv run python src/validate.py
"""

import numpy as np

from aerosim import naca4, Geometry, solve


def check(label, value, expected, tol):
    ok = abs(value - expected) <= tol
    flag = "OK " if ok else "FAIL"
    print(f"[{flag}] {label}: got {value:+.4f}, expected {expected:+.4f} (+/-{tol})")
    return ok


def main():
    passed = True

    # 1. Symmetric airfoil at zero incidence -> no lift.
    geom = Geometry(*naca4("0012", n_panels=200))
    s0 = solve(geom, 0.0)
    passed &= check("NACA0012 a=0 Cl", s0.cl, 0.0, 1e-3)

    # 2. Lift-curve slope. Thin-airfoil theory gives 2*pi/rad for a flat
    #    plate; finite thickness raises it by ~ (1 + 0.77 t/c).
    a1, a2 = 4.0, 8.0
    cl1 = solve(geom, a1).cl
    cl2 = solve(geom, a2).cl
    slope_per_deg = (cl2 - cl1) / (a2 - a1)
    slope_per_rad = slope_per_deg * 180.0 / np.pi
    expected_slope = 2 * np.pi * (1 + 0.77 * 0.12)
    passed &= check("NACA0012 dCl/da [1/rad]", slope_per_rad, expected_slope, 0.25)

    # 3. Cambered airfoil at a=0 has positive lift (zero-lift angle < 0).
    cam = Geometry(*naca4("2412", n_panels=200))
    sc = solve(cam, 0.0)
    passed &= check("NACA2412 a=0 Cl", sc.cl, 0.25, 0.06)

    # 3b. Cambered airfoil has a nose-down (negative) moment about c/4;
    #     a symmetric airfoil's quarter-chord moment is ~0.
    passed &= check("NACA2412 a=4 Cm_c/4", solve(cam, 4.0).cm_qc, -0.053, 0.02)
    passed &= check("NACA0012 a=5 Cm_c/4", solve(geom, 5.0).cm_qc, 0.0, 0.01)

    # 4. d'Alembert: pressure drag must vanish in potential flow.
    s5 = solve(geom, 5.0)
    passed &= check("NACA0012 a=5 Cd (should be ~0)", s5.cd, 0.0, 1e-2)

    # 5. Stagnation pressure: max Cp on the surface should reach ~1.
    passed &= check("NACA0012 a=5 max Cp", s5.cp.max(), 1.0, 0.03)

    # 6. Lift from circulation (Kutta-Joukowski) should match the pressure
    #    integral. Cl = -2 * Gamma / (Vinf * c); Gamma = gamma * perimeter
    #    (the minus sign reflects the CCW-positive circulation convention).
    cl_kj = -2.0 * s5.gamma * geom.perimeter / s5.vinf
    passed &= check("NACA0012 a=5 Cl (KJ vs pressure)", cl_kj, s5.cl, 0.02)

    print("\n" + ("All checks passed." if passed else "Some checks FAILED."))
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
