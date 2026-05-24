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
from fastapi import FastAPI, HTTPException
from fastapi.responses import FileResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel

from .airfoil import naca4
from .airfoil_io import list_samples, load_sample_text, parse_dat, repanel
from .flowfield import streamlines_from_grid, velocity_field
from .panel import Geometry, solve

STATIC_DIR = Path(__file__).parent / "static"

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


def _pack(geom: Geometry, name: str, alpha: float, re: float) -> dict:
    """Solve and assemble the JSON payload the front-end plots (source-agnostic)."""
    sol = solve(geom, alpha, re=re)
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
    half = geom.n // 2

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
        },
    }


class CustomRequest(BaseModel):
    """An uploaded .dat file (``dat``) or a bundled airfoil (``sample``)."""

    dat: str | None = None
    sample: str | None = None
    name: str | None = None
    alpha: float = 5.0
    re_log: float = 6.0
    panels: int = 160


def create_app() -> FastAPI:
    app = FastAPI(title="aerosim web", docs_url="/api/docs")

    @app.get("/api/solve")
    def api_solve(
        m: int = 2,
        p: int = 4,
        t: int = 12,
        alpha: float = 5.0,
        re_log: float = 6.0,
        panels: int = 160,
    ):
        """Solve a NACA 4-digit airfoil (from the sliders)."""
        code, p = _naca_code(m, p, t)
        geom = Geometry(*naca4(code, n_panels=panels))
        resp = _pack(geom, f"NACA {code}", alpha, 10.0**re_log)
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
            if req.sample:
                text = load_sample_text(req.sample)
            elif req.dat:
                text = req.dat
            else:
                raise HTTPException(400, "provide either 'dat' text or a 'sample' key")
            x, y, parsed_name = parse_dat(text)
            rx, ry = repanel(x, y, int(np.clip(req.panels, 40, 400)))
            geom = Geometry(rx, ry)
        except HTTPException:
            raise
        except Exception as exc:  # parse/geometry failure -> friendly 400
            raise HTTPException(400, f"could not load airfoil: {exc}")

        name = req.name or parsed_name or "loaded airfoil"
        resp = _pack(geom, name, req.alpha, 10.0**req.re_log)
        resp["source"] = "custom"
        resp["n_points"] = int(len(x))
        return resp

    @app.get("/")
    def index():
        return FileResponse(STATIC_DIR / "index.html")

    app.mount("/static", StaticFiles(directory=STATIC_DIR), name="static")
    return app


app = create_app()
