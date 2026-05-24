"""FastAPI backend for the web UI.

A thin JSON layer over the existing solver — no physics lives here. The browser
sends airfoil/flow parameters, this returns everything needed to draw the flow:
the airfoil outline, a pressure (``Cp``) field, streamline polylines, the
surface pressure distribution, and the force/drag coefficients (inviscid plus
the viscous boundary-layer estimate).

Run it via ``src/web.py`` (``uv run python src/web.py``).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
from fastapi import FastAPI
from fastapi.responses import FileResponse
from fastapi.staticfiles import StaticFiles

from .airfoil import naca4
from .flowfield import streamlines_from_grid, velocity_field
from .panel import Geometry, solve

STATIC_DIR = Path(__file__).parent / "static"

# Flow-field window and resolution (matches the desktop explorer's framing).
XLIM = (-0.6, 1.6)
YLIM = (-0.7, 0.7)
FIELD_NX, FIELD_NY = 130, 90


def _naca_code(m: int, p: int, t: int) -> tuple[str, int]:
    """Build a 4-digit code, mirroring the desktop UI's camber-position rule."""
    if m > 0 and p == 0:  # camber needs a non-zero position
        p = 1
    return f"{m}{p}{t:02d}", p


def _round_list(a, decimals=4):
    return np.round(np.asarray(a, dtype=float), decimals).tolist()


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
        """Solve at the given parameters and return everything the UI plots."""
        code, p = _naca_code(m, p, t)
        re = 10.0**re_log
        geom = Geometry(*naca4(code, n_panels=panels))
        sol = solve(geom, alpha, re=re)
        bl = sol.bl

        # Pressure field + streamlines on the background grid.
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
            "code": code,
            "p": int(p),
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

    @app.get("/")
    def index():
        return FileResponse(STATIC_DIR / "index.html")

    app.mount("/static", StaticFiles(directory=STATIC_DIR), name="static")
    return app


app = create_app()
