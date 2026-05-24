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
