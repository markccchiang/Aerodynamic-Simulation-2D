"""Hess-Smith constant-strength source + vortex panel method (2D, incompressible).

Each of the ``N`` panels carries a constant-strength source (one unknown per
panel) plus a single vortex strength shared by every panel. That gives ``N + 1``
unknowns, closed by ``N`` flow-tangency conditions (no flow through the surface)
and one Kutta condition (smooth flow off the trailing edge).

The induced-velocity formulas are the closed-form integrals of a point source
and a point vortex over a straight panel, evaluated in panel-local coordinates
and rotated back to global axes.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from matplotlib.path import Path

TWO_PI = 2.0 * np.pi


@dataclass
class Geometry:
    """Paneled airfoil geometry and precomputed panel quantities."""

    x: np.ndarray  # node x, shape (N+1,)
    y: np.ndarray  # node y, shape (N+1,)

    # populated in __post_init__
    n: int = field(init=False)
    xc: np.ndarray = field(init=False)
    yc: np.ndarray = field(init=False)

    def __post_init__(self):
        self.xa, self.ya = self.x[:-1], self.y[:-1]   # panel start nodes
        self.xb, self.yb = self.x[1:], self.y[1:]     # panel end nodes
        dx, dy = self.xb - self.xa, self.yb - self.ya
        self.length = np.hypot(dx, dy)
        self.theta = np.arctan2(dy, dx)
        self.ct, self.st = np.cos(self.theta), np.sin(self.theta)
        self.xc = 0.5 * (self.xa + self.xb)           # control points
        self.yc = 0.5 * (self.ya + self.yb)
        # Reference tangent (local +x) and normal (local +y).
        self.tx, self.ty = self.ct, self.st
        self.nx, self.ny = -self.st, self.ct
        self.n = len(self.length)
        self.perimeter = float(self.length.sum())
        self._exterior_side = self._detect_exterior_side()
        # Unit outward normal (points into the flow).
        self.noutx = self._exterior_side * self.nx
        self.nouty = self._exterior_side * self.ny

    def _detect_exterior_side(self) -> int:
        """Return +1 if the local +y normal points into the flow, else -1.

        Steps a small distance from each control point along the local +y
        normal and tests whether that point lies inside the airfoil polygon.
        A majority vote makes this robust near the thin trailing edge.
        """
        poly = Path(np.column_stack([self.x, self.y]))
        eps = 1e-4 * max(float(np.ptp(self.x)), 1e-9)
        probe = np.column_stack([self.xc + eps * self.nx, self.yc + eps * self.ny])
        inside = poly.contains_points(probe)
        # If +y normal lands inside, the flow side is -y -> exterior_side = -1.
        return -1 if inside.mean() > 0.5 else 1


def induced(px, py, geom: Geometry):
    """Global induced velocities at points ``(px, py)`` for unit singularities.

    Returns four ``(M, N)`` arrays ``us_x, us_y, uv_x, uv_y`` where ``M`` is the
    number of field points and ``N`` the number of panels. Column ``j`` holds the
    velocity induced by a unit-strength source (``us_*``) or vortex (``uv_*``) on
    panel ``j``; multiply by the panel strengths and sum to get the field.
    """
    px = np.atleast_1d(np.asarray(px, dtype=float))
    py = np.atleast_1d(np.asarray(py, dtype=float))

    dx = px[:, None] - geom.xa[None, :]
    dy = py[:, None] - geom.ya[None, :]
    ct, st = geom.ct[None, :], geom.st[None, :]

    # Field point in panel-local coordinates (panel along local +x from 0..L).
    xl = dx * ct + dy * st
    yl = -dx * st + dy * ct
    L = geom.length[None, :]

    r1sq = xl**2 + yl**2
    r2sq = (xl - L) ** 2 + yl**2
    lnr = 0.5 * np.log(r1sq / r2sq)                      # ln(r1 / r2)
    dth = np.arctan2(yl, xl - L) - np.arctan2(yl, xl)    # theta2 - theta1

    # Local-frame velocity components (unit strength).
    up_s, wp_s = lnr / TWO_PI, dth / TWO_PI              # source
    up_v, wp_v = -dth / TWO_PI, lnr / TWO_PI             # vortex (source rot. 90 deg)

    # Rotate local -> global.
    us_x = up_s * ct - wp_s * st
    us_y = up_s * st + wp_s * ct
    uv_x = up_v * ct - wp_v * st
    uv_y = up_v * st + wp_v * ct
    return us_x, us_y, uv_x, uv_y


@dataclass
class Solution:
    geom: Geometry
    alpha: float            # angle of attack, degrees
    vinf: float
    sigma: np.ndarray       # source strengths, (N,)
    gamma: float            # shared vortex strength
    cp: np.ndarray          # surface pressure coefficient at control points
    vt: np.ndarray          # surface tangential velocity at control points
    cl: float               # lift coefficient
    cd: float               # pressure drag (≈0; a numerical accuracy check)
    cm_qc: float            # moment coefficient about the quarter chord


def solve(geom: Geometry, alpha_deg: float, vinf: float = 1.0) -> Solution:
    """Solve the panel system at the given angle of attack."""
    a = np.radians(alpha_deg)
    uinf, winf = vinf * np.cos(a), vinf * np.sin(a)
    n = geom.n

    us_x, us_y, uv_x, uv_y = induced(geom.xc, geom.yc, geom)

    # Self-influence: the limit approaching each control point from the flow
    # side. Source -> 0.5 outward-normal velocity; vortex -> 0.5 tangential.
    s = geom._exterior_side
    di = np.arange(n)
    us_x[di, di] = -(0.5 * s) * geom.st
    us_y[di, di] = (0.5 * s) * geom.ct
    uv_x[di, di] = -(0.5 * s) * geom.ct
    uv_y[di, di] = -(0.5 * s) * geom.st

    nx, ny = geom.nx, geom.ny
    tx, ty = geom.tx, geom.ty

    A = np.zeros((n + 1, n + 1))
    b = np.zeros(n + 1)

    # Flow tangency (no penetration) at every control point.
    A[:n, :n] = us_x * nx[:, None] + us_y * ny[:, None]
    A[:n, n] = (uv_x * nx[:, None] + uv_y * ny[:, None]).sum(axis=1)
    b[:n] = -(uinf * nx + winf * ny)

    # Kutta condition: equal-and-opposite tangential velocity on the two
    # trailing-edge panels (the first and last in the node ordering).
    i0, i1 = 0, n - 1
    A[n, :n] = (us_x[i0] * tx[i0] + us_y[i0] * ty[i0]) + (
        us_x[i1] * tx[i1] + us_y[i1] * ty[i1]
    )
    A[n, n] = (uv_x[i0] * tx[i0] + uv_y[i0] * ty[i0]).sum() + (
        uv_x[i1] * tx[i1] + uv_y[i1] * ty[i1]
    ).sum()
    b[n] = -(
        (uinf * tx[i0] + winf * ty[i0]) + (uinf * tx[i1] + winf * ty[i1])
    )

    sol = np.linalg.solve(A, b)
    sigma, gamma = sol[:n], float(sol[n])

    # Surface velocity and pressure at the control points.
    vx = us_x @ sigma + gamma * uv_x.sum(axis=1) + uinf
    vy = us_y @ sigma + gamma * uv_y.sum(axis=1) + winf
    vt = vx * tx + vy * ty
    cp = 1.0 - (vt / vinf) ** 2

    # Integrate pressure over the surface (force = -∮ Cp n_out ds, chord = 1).
    fx = -np.sum(cp * geom.noutx * geom.length)
    fy = -np.sum(cp * geom.nouty * geom.length)
    cl = -fx * np.sin(a) + fy * np.cos(a)
    cd = fx * np.cos(a) + fy * np.sin(a)          # ~0 by d'Alembert's paradox
    # Aerodynamic moment convention: positive = nose-up (= -CCW moment).
    cm_qc = np.sum(
        cp * geom.length * ((geom.xc - 0.25) * geom.nouty - geom.yc * geom.noutx)
    )

    return Solution(geom, alpha_deg, vinf, sigma, gamma, cp, vt, cl, cd, cm_qc)
