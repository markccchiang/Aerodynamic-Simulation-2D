"""Reconstruct the velocity field around a solved airfoil for streamlines."""

from __future__ import annotations

import numpy as np
from matplotlib.path import Path

from .panel import Solution, induced


def velocity_field(sol: Solution, xs: np.ndarray, ys: np.ndarray, mask_body: bool = True):
    """Evaluate the flow velocity on a grid spanned by ``xs`` x ``ys``.

    Returns ``(X, Y, U, V)`` meshgrids. Points inside the airfoil are set to NaN
    when ``mask_body`` is True so they are skipped by ``streamplot``.
    """
    X, Y = np.meshgrid(xs, ys)
    px, py = X.ravel(), Y.ravel()

    us_x, us_y, uv_x, uv_y = induced(px, py, sol.geom)
    a = np.radians(sol.alpha)
    U = us_x @ sol.sigma + sol.gamma * uv_x.sum(axis=1) + sol.vinf * np.cos(a)
    V = us_y @ sol.sigma + sol.gamma * uv_y.sum(axis=1) + sol.vinf * np.sin(a)
    U = U.reshape(X.shape)
    V = V.reshape(X.shape)

    if mask_body:
        poly = Path(np.column_stack([sol.geom.x, sol.geom.y]))
        inside = poly.contains_points(np.column_stack([px, py])).reshape(X.shape)
        U[inside] = np.nan
        V[inside] = np.nan

    return X, Y, U, V


def speed_and_cp(U: np.ndarray, V: np.ndarray, vinf: float = 1.0):
    """Return flow speed and pressure coefficient fields from velocity grids."""
    speed = np.hypot(U, V)
    cp = 1.0 - (speed / vinf) ** 2
    return speed, cp


def streamlines_from_grid(xs, ys, U, V, n_lines=46, ds=0.02, max_steps=800):
    """Integrate streamlines across a velocity grid (bilinear interp + RK2).

    Seeds ``n_lines`` points evenly along the inlet (left edge) and marches
    each one downstream with a constant arc-length step ``ds`` (so the polyline
    vertices are evenly spaced, independent of local speed). A line stops when
    it leaves the grid or hits a NaN cell — i.e. the masked airfoil body.

    Returns a list of polylines, each a dict with ``x`` and ``y`` lists. This
    is the matplotlib-free equivalent of ``streamplot`` for the web front-end.
    """
    xs = np.asarray(xs, dtype=float)
    ys = np.asarray(ys, dtype=float)
    nx, ny = len(xs), len(ys)
    x0, x1, y0, y1 = xs[0], xs[-1], ys[0], ys[-1]
    dx = (x1 - x0) / (nx - 1)
    dy = (y1 - y0) / (ny - 1)

    def sample(px, py):
        """Bilinearly interpolate (u, v) at a point; NaN outside / in body."""
        fx, fy = (px - x0) / dx, (py - y0) / dy
        i, j = int(np.floor(fx)), int(np.floor(fy))
        if i < 0 or i >= nx - 1 or j < 0 or j >= ny - 1:
            return np.nan, np.nan
        tx, ty = fx - i, fy - j

        def bil(F):
            return (
                (F[j, i] * (1 - tx) + F[j, i + 1] * tx) * (1 - ty)
                + (F[j + 1, i] * (1 - tx) + F[j + 1, i + 1] * tx) * ty
            )

        return bil(U), bil(V)

    margin = 0.02 * (y1 - y0)
    lines = []
    for sy in np.linspace(y0 + margin, y1 - margin, n_lines):
        px, py = x0 + 1e-6, float(sy)
        xl, yl = [px], [py]
        for _ in range(max_steps):
            u1, v1 = sample(px, py)
            sp1 = np.hypot(u1, v1)
            if not np.isfinite(sp1) or sp1 < 1e-6:
                break
            # RK2 (midpoint), stepping by unit-speed direction.
            mx, my = px + 0.5 * ds * u1 / sp1, py + 0.5 * ds * v1 / sp1
            u2, v2 = sample(mx, my)
            sp2 = np.hypot(u2, v2)
            if not np.isfinite(sp2) or sp2 < 1e-6:
                break
            px, py = px + ds * u2 / sp2, py + ds * v2 / sp2
            if not (x0 <= px <= x1 and y0 <= py <= y1):
                break
            xl.append(px)
            yl.append(py)
        if len(xl) > 3:
            lines.append({"x": xl, "y": yl})
    return lines
