"""FastAPI backend for the web UI.

A thin JSON layer over the existing solver — no physics lives here. The browser
sends airfoil/flow parameters, this returns everything needed to draw the flow:
the airfoil outline, a pressure (``Cp``) field, streamline polylines, the
surface pressure distribution, and the force/drag coefficients (inviscid plus
the viscous boundary-layer estimate).

Airfoils come from three sources, all packed into the same response shape:
``GET /api/solve`` (NACA 4-digit from sliders), ``GET /api/samples`` +
``POST /api/solve_custom`` with a ``sample`` key (bundled ``.dat``), or
``POST /api/solve_custom`` with ``dat`` text (an uploaded coordinate file).

Run it via ``src/web.py`` (``uv run python src/web.py``).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
from fastapi import FastAPI, HTTPException, Query
from fastapi.responses import FileResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel, Field

from .airfoil import naca4
from .airfoil_io import list_samples, load_sample_text, parse_dat, repanel
from .flowfield import streamlines_from_grid, velocity_field
from .panel import Geometry, solve

STATIC_DIR = Path(__file__).parent / "static"

# Request bounds, shared by the GET query params and the POST request models so
# every endpoint rejects the same things. Out-of-range input is a 422 from the
# validation layer rather than a 500 from deep inside the solver; ``panels`` is
# capped because the panel matrix is (N+1)^2 and the server would otherwise
# happily try to allocate gigabytes for it.
ALPHA_MIN, ALPHA_MAX = -30.0, 30.0
RE_LOG_MIN, RE_LOG_MAX = 4.0, 9.0
PANELS_MIN, PANELS_MAX = 40, 400

# Flow-field window and resolution for the streamline/Cp background grid.
XLIM = (-0.6, 1.6)
YLIM = (-0.7, 0.7)
FIELD_NX, FIELD_NY = 130, 90


def _naca_code(m: int, p: int, t: int) -> tuple[str, int]:
    """Build a 4-digit code; camber needs a non-zero position, so snap P 0->1."""
    if m > 0 and p == 0:  # camber needs a non-zero position
        p = 1
    return f"{m}{p}{t:02d}", p


def _round_list(a, decimals=4):
    return np.round(np.asarray(a, dtype=float), decimals).tolist()


def _pack(geom: Geometry, name: str, alpha: float, re: float,
          couple: bool = False) -> dict:
    """Solve and assemble the JSON payload the front-end plots (source-agnostic).

    ``couple=True`` runs the two-way viscous-inviscid coupling, so ``cl``/``cp``
    react to viscosity; otherwise the boundary-layer estimate is uncoupled.
    """
    sol = solve(geom, alpha, re=re, couple=couple)
    bl = sol.bl

    xs = np.linspace(*XLIM, FIELD_NX)
    ys = np.linspace(*YLIM, FIELD_NY)
    _, _, U, V = velocity_field(sol, xs, ys)
    cp_field = 1.0 - (np.hypot(U, V) / sol.vinf) ** 2
    lines = streamlines_from_grid(xs, ys, U, V)

    # JSON has no NaN, so masked (in-body) cells become null.
    cp_rows = [
        [None if not np.isfinite(v) else round(float(v), 3) for v in row]
        for row in cp_field
    ]
    # Split the surface runs at the actual leading-edge control point. The LE
    # is only at the mid index for a symmetric section; elsewhere it sits 1-2
    # panels off, which mis-assigns Cp points right at the suction peak.
    half = int(np.argmin(geom.xc)) + 1

    def opt(x):  # finite float or None (transition may be absent)
        return None if not np.isfinite(x) else round(float(x), 3)

    return {
        "name": name,
        "alpha": float(alpha),
        "re": float(re),
        "geom": {"x": _round_list(geom.x, 5), "y": _round_list(geom.y, 5)},
        "field": {"x": _round_list(xs, 4), "y": _round_list(ys, 4), "cp": cp_rows},
        "streamlines": [
            {"x": _round_list(ln["x"], 4), "y": _round_list(ln["y"], 4)}
            for ln in lines
        ],
        "surface": {
            "x_upper": _round_list(geom.xc[:half], 4),
            "cp_upper": _round_list(sol.cp[:half], 4),
            "x_lower": _round_list(geom.xc[half:], 4),
            "cp_lower": _round_list(sol.cp[half:], 4),
        },
        "coeffs": {
            "cl": round(float(sol.cl), 4),
            "cd_pressure": round(float(sol.cd), 5),
            "cd_visc": round(float(bl.cd), 5),
            "cm": round(float(sol.cm_qc), 4),
            "x_tr_upper": opt(bl.upper.x_transition),
            "x_tr_lower": opt(bl.lower.x_transition),
            "separated": bool(bl.separated),
            "coupled": bool(sol.coupled),
            "converged": bool(sol.converged),
            "n_iter": int(sol.n_iter),
        },
    }


def _custom_geom(dat, sample, panels):
    """Build a Geometry from uploaded ``dat`` text or a bundled ``sample`` key.

    Returns ``(geom, n_points, parsed_name)``. Raises ``HTTPException(400)`` if
    neither input is given; let parse/geometry errors propagate to the caller.
    """
    if sample:
        text = load_sample_text(sample)
    elif dat:
        text = dat
    else:
        raise HTTPException(400, "provide either 'dat' text or a 'sample' key")
    x, y, parsed_name = parse_dat(text)
    rx, ry = repanel(x, y, int(np.clip(panels, 40, 400)))
    return Geometry(rx, ry), int(len(x)), parsed_name


def _polar(geom: Geometry, name: str, re: float, couple: bool,
           amin: float, amax: float, astep: float) -> dict:
    """Sweep angle of attack into lift-curve / drag-polar arrays.

    Always returns the ``uncoupled`` branch (inviscid ``cl`` with the uncoupled
    profile ``cd``); when ``couple`` is set it also returns a ``coupled`` branch
    so the front-end can overlay the viscous decambering. The step is floored and
    the point count capped to keep a stray request from launching a huge sweep.
    """
    astep = max(abs(astep), 0.25)
    n = min(max(int(round((amax - amin) / astep)) + 1, 2), 200)
    alphas = [round(amin + i * astep, 3) for i in range(n)]
    base = {"cl": [], "cd": [], "cm": []}
    cpl = {"cl": [], "cd": [], "cm": [], "converged": []}
    for a in alphas:
        s = solve(geom, a, re=re)
        base["cl"].append(round(float(s.cl), 4))
        base["cd"].append(round(float(s.cd_visc), 5))
        base["cm"].append(round(float(s.cm_qc), 4))
        if couple:
            c = solve(geom, a, re=re, couple=True)
            cpl["cl"].append(round(float(c.cl), 4))
            cpl["cd"].append(round(float(c.cd_visc), 5))
            cpl["cm"].append(round(float(c.cm_qc), 4))
            cpl["converged"].append(bool(c.converged))
    return {
        "name": name,
        "re": float(re),
        "alpha": alphas,
        "uncoupled": base,
        "coupled": cpl if couple else None,
    }


class CustomRequest(BaseModel):
    """An uploaded .dat file (``dat``) or a bundled airfoil (``sample``)."""

    dat: str | None = None
    sample: str | None = None
    name: str | None = None
    alpha: float = Field(5.0, ge=ALPHA_MIN, le=ALPHA_MAX)
    re_log: float = Field(6.0, ge=RE_LOG_MIN, le=RE_LOG_MAX)
    panels: int = Field(160, ge=PANELS_MIN, le=PANELS_MAX)
    couple: bool = False


class PolarRequest(BaseModel):
    """Sweep request for a loaded airfoil's lift curve / drag polar."""

    dat: str | None = None
    sample: str | None = None
    name: str | None = None
    re_log: float = Field(6.0, ge=RE_LOG_MIN, le=RE_LOG_MAX)
    panels: int = Field(160, ge=PANELS_MIN, le=PANELS_MAX)
    couple: bool = False
    amin: float = Field(-6.0, ge=ALPHA_MIN, le=ALPHA_MAX)
    amax: float = Field(14.0, ge=ALPHA_MIN, le=ALPHA_MAX)
    astep: float = Field(1.0, ge=0.25, le=10.0)


def create_app() -> FastAPI:
    app = FastAPI(title="aerosim web", docs_url="/api/docs")

    @app.get("/api/solve")
    def api_solve(
        m: int = Query(2, ge=0, le=9),
        p: int = Query(4, ge=0, le=9),
        t: int = Query(12, ge=1, le=40),
        alpha: float = Query(5.0, ge=ALPHA_MIN, le=ALPHA_MAX),
        re_log: float = Query(6.0, ge=RE_LOG_MIN, le=RE_LOG_MAX),
        panels: int = Query(160, ge=PANELS_MIN, le=PANELS_MAX),
        couple: bool = False,
    ):
        """Solve a NACA 4-digit airfoil (from the sliders)."""
        code, p = _naca_code(m, p, t)
        geom = Geometry(*naca4(code, n_panels=panels))
        resp = _pack(geom, f"NACA {code}", alpha, 10.0**re_log, couple)
        resp["source"] = "naca"
        resp["code"] = code
        resp["p"] = int(p)
        return resp

    @app.get("/api/samples")
    def api_samples():
        """List the bundled sample airfoils for the dropdown."""
        return {"samples": [{"key": k, "name": n} for k, n in list_samples()]}

    @app.post("/api/solve_custom")
    def api_solve_custom(req: CustomRequest):
        """Solve a loaded airfoil: a bundled ``sample`` or uploaded ``dat`` text."""
        try:
            geom, n_points, parsed_name = _custom_geom(req.dat, req.sample, req.panels)
        except HTTPException:
            raise
        except Exception as exc:  # parse/geometry failure -> friendly 400
            raise HTTPException(400, f"could not load airfoil: {exc}")

        name = req.name or parsed_name or "loaded airfoil"
        resp = _pack(geom, name, req.alpha, 10.0**req.re_log, req.couple)
        resp["source"] = "custom"
        resp["n_points"] = n_points
        return resp

    @app.get("/api/polar")
    def api_polar(
        m: int = Query(2, ge=0, le=9),
        p: int = Query(4, ge=0, le=9),
        t: int = Query(12, ge=1, le=40),
        re_log: float = Query(6.0, ge=RE_LOG_MIN, le=RE_LOG_MAX),
        panels: int = Query(160, ge=PANELS_MIN, le=PANELS_MAX),
        couple: bool = False,
        amin: float = Query(-6.0, ge=ALPHA_MIN, le=ALPHA_MAX),
        amax: float = Query(14.0, ge=ALPHA_MIN, le=ALPHA_MAX),
        astep: float = Query(1.0, ge=0.25, le=10.0),
    ):
        """Lift curve / drag polar for a NACA 4-digit airfoil (alpha sweep)."""
        code, _p = _naca_code(m, p, t)
        geom = Geometry(*naca4(code, n_panels=panels))
        return _polar(geom, f"NACA {code}", 10.0**re_log, couple, amin, amax, astep)

    @app.post("/api/polar_custom")
    def api_polar_custom(req: PolarRequest):
        """Lift curve / drag polar for a loaded airfoil (alpha sweep)."""
        try:
            geom, _n, parsed_name = _custom_geom(req.dat, req.sample, req.panels)
        except HTTPException:
            raise
        except Exception as exc:  # parse/geometry failure -> friendly 400
            raise HTTPException(400, f"could not load airfoil: {exc}")

        name = req.name or parsed_name or "loaded airfoil"
        return _polar(geom, name, 10.0**req.re_log, req.couple, req.amin, req.amax, req.astep)

    @app.get("/")
    def index():
        return FileResponse(STATIC_DIR / "index.html")

    app.mount("/static", StaticFiles(directory=STATIC_DIR), name="static")
    return app


app = create_app()
