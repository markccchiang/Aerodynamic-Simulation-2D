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


def check_true(label, cond, detail=""):
    flag = "OK " if cond else "FAIL"
    print(f"[{flag}] {label}{(': ' + detail) if detail else ''}")
    return bool(cond)


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

    # 7. Exact reference case. Potential flow past a circle has the closed-form
    #    solution Cp = 1 - 4 sin^2(theta). Unlike every other check in this file
    #    that is analysis rather than a trend or a published number, so it pins
    #    induced() and the panel assembly directly. Built in code at full
    #    precision: the bundled circle.dat stores 6 decimals, which limits the
    #    same comparison to ~4e-4 (see check 19).
    nc = 480
    tc = np.linspace(0.0, 2.0 * np.pi, nc + 1)
    cx, cy = 0.5 + 0.5 * np.cos(tc), 0.5 * np.sin(tc)
    cx[0] = cx[-1] = 1.0
    cy[0] = cy[-1] = 0.0
    cy[nc // 2] = 0.0
    cyl = Geometry(cx, cy)
    cyl_cp = solve(cyl, 0.0).cp
    cyl_th = np.mod(np.arctan2(cyl.yc, cyl.xc - 0.5), 2.0 * np.pi)
    cyl_err = float(np.max(np.abs(cyl_cp - (1.0 - 4.0 * np.sin(cyl_th) ** 2))))
    passed &= check_true("Circle Cp == 1 - 4 sin^2(theta) (exact)", cyl_err < 1e-10,
                         f"max |dCp| = {cyl_err:.2e}")

    # ---- Viscous boundary-layer correction (uncoupled profile drag) ----
    print("\n-- viscous boundary-layer correction --")

    # 8. Profile drag of NACA0012 at a=0. Published Cd0 (free transition) is
    #    ~0.0080 at Re=1e6; this integral method should land close.
    v10 = solve(geom, 0.0, re=1e6)
    passed &= check("NACA0012 a=0 Re=1e6 Cd (profile)", v10.bl.cd, 0.0080, 0.0025)

    # 9. The viscous Cd is real (positive) drag -- the whole point of the
    #    correction. The inviscid pressure Cd stays ~0 (d'Alembert) regardless.
    passed &= check_true("Viscous Cd is positive", v10.bl.cd > 0,
                         f"Cd={v10.bl.cd:.4f}")

    # 10. Drag falls as Reynolds number rises (thinner boundary layer).
    v90 = solve(geom, 0.0, re=9e6)
    passed &= check_true("Cd decreases with Reynolds number", v90.bl.cd < v10.bl.cd,
                         f"Cd(9e6)={v90.bl.cd:.4f} < Cd(1e6)={v10.bl.cd:.4f}")

    # 11. Drag rises with incidence (more suction-side adverse gradient).
    v18 = solve(geom, 8.0, re=1e6)
    passed &= check_true("Cd increases with angle of attack", v18.bl.cd > v10.bl.cd,
                         f"Cd(a=8)={v18.bl.cd:.4f} > Cd(a=0)={v10.bl.cd:.4f}")

    # 12. Transition moves toward the leading edge as Reynolds number rises.
    passed &= check_true("Transition advances with Reynolds number",
                         v90.bl.upper.x_transition < v10.bl.upper.x_transition,
                         f"x_tr(9e6)={v90.bl.upper.x_transition:.3f} < "
                         f"x_tr(1e6)={v10.bl.upper.x_transition:.3f}")

    # 13. Integrated skin friction is a positive subset of the profile drag.
    passed &= check_true("0 < friction drag < profile drag",
                         0.0 < v10.bl.cd_friction < v10.bl.cd,
                         f"Cdf={v10.bl.cd_friction:.4f}, Cd={v10.bl.cd:.4f}")

    # 14. Attaching a viscous estimate must not perturb the inviscid solution.
    passed &= check("Inviscid Cl unchanged by re=", v10.cl, s0.cl, 1e-12)
    passed &= check("Inviscid Cd unchanged by re=", v10.cd, s0.cd, 1e-12)

    # 15. Trailing-edge separation onset: attached at a=0, separated near stall.
    v_stall = solve(geom, 14.0, re=1e6)
    passed &= check_true("Attached at a=0, separates near stall (a=14)",
                         (not v10.bl.separated) and v_stall.bl.separated,
                         f"x_sep(a=14)={v_stall.bl.upper.x_separation:.3f}")

    # ---- Viscous-inviscid coupling (two-way; couple=True) ----
    print("\n-- viscous-inviscid coupling --")

    # C1. Coupling needs a Reynolds number to march the boundary layer.
    try:
        solve(geom, 4.0, couple=True)
        raised = False
    except ValueError:
        raised = True
    passed &= check_true("couple=True without re= raises", raised)

    # C2. Two-way coupling reduces lift (viscous decambering) versus the inviscid
    #     solve, and the iteration converges for an attached case. The uncoupled
    #     default (couple=False) still returns the inviscid cl unchanged.
    a_c = 5.0
    inv = solve(cam, a_c, re=1e6)                 # uncoupled: cl == inviscid
    cpl = solve(cam, a_c, re=1e6, couple=True)    # two-way coupled
    passed &= check_true("Coupled solve converged", cpl.converged,
                         f"n_iter={cpl.n_iter}")
    passed &= check_true("Coupling reduces lift (decambering)",
                         0.02 < (inv.cl - cpl.cl) < 0.20,
                         f"cl {inv.cl:.4f} -> {cpl.cl:.4f} (dCl={cpl.cl - inv.cl:+.4f})")

    # C3. The suction peak softens under coupling (least-negative Cp rises).
    passed &= check_true("Coupling softens the suction peak",
                         cpl.cp.min() > inv.cp.min(),
                         f"min Cp {inv.cp.min():.3f} -> {cpl.cp.min():.3f}")

    # C4. A symmetric airfoil at zero lift gains no spurious lift from coupling.
    sym_c = solve(geom, 0.0, re=1e6, couple=True)
    passed &= check("Coupled symmetric a=0 Cl ~ 0", sym_c.cl, 0.0, 1e-3)

    # C5. The lift loss grows as Reynolds number falls (thicker boundary layer).
    d_hi = solve(cam, 4.0, re=3e6).cl - solve(cam, 4.0, re=3e6, couple=True).cl
    d_lo = solve(cam, 4.0, re=1e6).cl - solve(cam, 4.0, re=1e6, couple=True).cl
    passed &= check_true("Lift loss grows as Re falls", d_lo > d_hi > 0,
                         f"dCl(1e6)={d_lo:.4f} > dCl(3e6)={d_hi:.4f}")

    # C6. Profile drag is still a positive estimate after coupling.
    passed &= check_true("Coupled profile Cd positive", cpl.cd_visc > 0,
                         f"cd_visc={cpl.cd_visc:.4f}")

    # ---- Loading airfoils from .dat coordinate files ----
    print("\n-- .dat airfoil loading --")
    from aerosim.airfoil_io import (
        parse_dat, repanel, format_selig, format_lednicer, list_samples, load_sample_text,
    )

    x2, y2 = naca4("2412", n_panels=200)
    cl_ref = solve(Geometry(x2, y2), 4.0).cl  # direct NACA, the reference

    # 16. Selig round-trip: write coords -> parse -> re-panel -> solve == direct.
    sx, sy, sname = parse_dat(format_selig(x2, y2, "NACA 2412"))
    passed &= check("Selig .dat round-trip Cl", solve(Geometry(*repanel(sx, sy, 200)), 4.0).cl,
                    cl_ref, 0.01)
    passed &= check_true("Selig title line parsed", sname == "NACA 2412", repr(sname))

    # 17. Lednicer format parses to the same geometry.
    lx, ly, _ = parse_dat(format_lednicer(x2, y2, "NACA 2412"))
    passed &= check("Lednicer .dat round-trip Cl", solve(Geometry(*repanel(lx, ly, 200)), 4.0).cl,
                    cl_ref, 0.01)

    # 18. Coordinates are normalized to unit chord (scale + offset invariant).
    nx, ny, _ = parse_dat(format_selig(x2 * 4.0 + 3.0, y2 * 4.0, "scaled"))
    passed &= check("Scaled/offset .dat normalized Cl", solve(Geometry(*repanel(nx, ny, 200)), 4.0).cl,
                    cl_ref, 0.01)

    # 19. The bundled circle.dat is that same analytic case on disk, so reading
    #     it back must still match Cp = 1 - 4 sin^2(theta). The tolerance is set
    #     by the file's 6-decimal coordinates, not by the solver.
    dx, dy, _ = parse_dat(load_sample_text("circle"))
    dg = Geometry(dx, dy)
    d_th = np.mod(np.arctan2(dg.yc, dg.xc - 0.5), 2.0 * np.pi)
    d_err = float(np.max(np.abs(solve(dg, 0.0).cp - (1.0 - 4.0 * np.sin(d_th) ** 2))))
    passed &= check_true("circle.dat Cp matches the analytic circle", d_err < 2e-3,
                         f"max |dCp| = {d_err:.2e}")

    # 20. Every bundled sample loads, re-panels, and solves sanely.
    samples = list_samples()
    passed &= check_true("Bundled samples present", len(samples) >= 1,
                         str([k for k, _ in samples]))
    all_ok = True
    for key, _name in samples:
        gx, gy = repanel(*parse_dat(load_sample_text(key))[:2], 200)
        s = solve(Geometry(gx, gy), 5.0, re=1e6)
        all_ok &= bool(np.isfinite(s.cl) and s.bl.cd > 0)
    passed &= check_true("All bundled samples solve (finite Cl, Cd>0)", all_ok)

    # ---- The layers above the physics: payload shape, request validation ----
    print("\n-- web payload & solver plumbing --")
    from pydantic import ValidationError

    from aerosim.flowfield import velocity_field
    from aerosim.webapp import CustomRequest, _naca_code, _pack

    # 21. The Cp payload must split the surface at the true leading edge. The LE
    #     only lands on the mid index for a symmetric section, so splitting at
    #     n//2 puts lower-surface points on the upper curve at the suction peak.
    split_ok = True
    for key, _name in samples:
        g = Geometry(*repanel(*parse_dat(load_sample_text(key))[:2], 160))
        surf = _pack(g, key, 5.0, 1e6)["surface"]
        le = int(np.argmin(g.xc))
        split_ok &= (
            len(surf["x_upper"]) == le + 1
            and len(surf["x_lower"]) == g.n - le - 1
            and np.allclose(surf["x_upper"], np.round(g.xc[: le + 1], 4))
            and np.allclose(surf["x_lower"], np.round(g.xc[le + 1 :], 4))
        )
    passed &= check_true("Cp payload splits at the leading edge", split_ok)

    # 22. Out-of-range requests are rejected by the validation layer instead of
    #     reaching the solver and surfacing as a 500.
    bad = (dict(alpha=999.0), dict(re_log=400.0), dict(re_log=-400.0),
           dict(panels=0), dict(panels=20000))
    rejected = 0
    for kw in bad:
        try:
            CustomRequest(sample="naca0012", **kw)
        except ValidationError:
            rejected += 1
    passed &= check_true("Out-of-range requests rejected", rejected == len(bad),
                         f"{rejected}/{len(bad)}")

    # 23. Camber needs a non-zero position, so the API snaps P 0 -> 1; without it
    #     M>0 with P=0 silently returns a symmetric section.
    code_snapped, p_snapped = _naca_code(2, 0, 12)
    passed &= check_true("Cambered NACA code snaps P 0 -> 1",
                         code_snapped == "2112" and p_snapped == 1, code_snapped)

    # 24. The panel factorization is cached on the Geometry, so a reused geometry
    #     must give bit-identical results to a freshly built one at every alpha.
    gx, gy = naca4("2412", n_panels=160)
    shared = Geometry(gx, gy)
    passed &= check_true(
        "Cached geometry matches a fresh solve",
        all(solve(shared, a).cl == solve(Geometry(gx, gy), a).cl
            for a in (-4.0, 0.0, 6.0, 11.0)),
    )

    # 25. The reconstructed flow field relaxes to the free stream far upstream
    #     and masks points inside the body (the web UI relies on both).
    ff = solve(Geometry(gx, gy), 6.0)
    _, _, U, V = velocity_field(ff, np.array([-30.0, 0.5]), np.array([0.0]))
    a_ff = np.radians(6.0)
    passed &= check_true(
        "Far field -> free stream, body masked",
        abs(U[0, 0] - np.cos(a_ff)) < 5e-3 and abs(V[0, 0] - np.sin(a_ff)) < 5e-3
        and bool(np.isnan(U[0, 1])) and bool(np.isnan(V[0, 1])),
        f"far (u,v)=({U[0, 0]:.4f}, {V[0, 0]:.4f}), in-body NaN={np.isnan(U[0, 1])}",
    )

    print("\n" + ("All checks passed." if passed else "Some checks FAILED."))
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
