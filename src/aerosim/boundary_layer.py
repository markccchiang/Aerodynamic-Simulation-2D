"""Viscous boundary-layer correction (uncoupled by default; optional coupling).

Given an inviscid :class:`~aerosim.panel.Solution`, march an *integral*
boundary-layer method along each surface from the stagnation point to the
trailing edge and estimate the **profile (viscous) drag** the potential-flow
model cannot produce. On its own this is a one-way correction: it reads the
inviscid edge velocities but does **not** feed back into the panel solve, so
lift, ``Cp`` and the d'Alembert ``cd`` check are left exactly as they were.

For *two-way* viscous-inviscid coupling, :func:`transpiration_velocity` turns
the displacement thickness into an equivalent wall-blowing velocity that
``panel.solve(..., couple=True)`` feeds back into the flow-tangency boundary
condition and iterates to convergence. Coupling then *does* change ``cl``/``cp``
(viscous decambering); see :mod:`aerosim.panel`.

Method (the classic "XFOIL-lite" chain):

* **Laminar run**   – Thwaites' integral method.
* **Transition**    – Michel's criterion (or forced at laminar separation).
* **Turbulent run** – Head's entrainment method with Ludwieg–Tillmann ``Cf``.
* **Profile drag**  – the Squire–Young formula evaluated at each trailing edge,
  which extrapolates the trailing-edge momentum thickness through the wake.

Everything is non-dimensionalised by chord (``c = 1``) and free-stream speed
(``Vinf``), so the kinematic viscosity is simply ``nu = 1 / Re`` and the edge
velocity is ``Ue / Vinf``. Reynolds number is chord-based: ``Re = Vinf c / nu``.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:  # avoid an import cycle (panel imports this module)
    from .panel import Solution

# Thwaites' laminar separation parameter.
_THWAITES_SEP = -0.09
# Head's method: entrainment shape factor H1 falls toward ~3.3 at separation;
# H1 ~ 3.5 corresponds to H ~ 2.6, which we treat as turbulent separation.
_H1_SEP = 3.5
_TURB_SEP_H = 2.6
# Separation within the last few percent of chord is a cusped-trailing-edge
# artifact of the integral method, not a real separated region; only count
# separation forward of this station as a genuine "separated" flag.
_SEP_X_LIMIT = 0.90
# Coupling: the same cusped-TE artifact makes delta*/Ue (and hence the
# transpiration velocity) spike over the last few percent chord. Taper the
# coupled wall-blowing smoothly to zero from this station to the TE so that
# artifact is not fed back into the inviscid solve (it diverges the iteration).
_TE_TAPER_X = 0.95


# --------------------------------------------------------------------- results


@dataclass
class SurfaceBL:
    """Boundary-layer state along one surface (stagnation -> trailing edge)."""

    name: str               # "upper" or "lower"
    s: np.ndarray           # arc length from the stagnation point
    x: np.ndarray           # control-point x/c (for plotting against the body)
    ue: np.ndarray          # edge velocity Ue / Vinf
    theta: np.ndarray       # momentum thickness theta / c
    delta_star: np.ndarray  # displacement thickness delta* / c
    H: np.ndarray           # shape factor delta* / theta
    cf: np.ndarray          # local skin-friction coefficient (ref. to Ue)
    i_transition: int       # index of the transition point (-1 if fully laminar)
    x_transition: float     # transition location x/c (nan if fully laminar)
    x_separation: float     # turbulent separation location x/c (nan if attached)
    separated: bool         # True if separation occurs forward of the TE cusp
    cd: float               # Squire–Young profile-drag contribution


@dataclass
class BoundaryLayer:
    """Combined viscous estimate from both surfaces."""

    re: float
    upper: SurfaceBL
    lower: SurfaceBL
    cd: float               # total profile drag (upper.cd + lower.cd)
    cd_friction: float      # integrated skin-friction drag (a diagnostic subset)
    separated: bool = field(init=False)

    def __post_init__(self):
        self.separated = self.upper.separated or self.lower.separated


# ------------------------------------------------------------- correlations


def _thwaites_l_H(lam: float):
    """Thwaites shear ``l(lambda)`` and shape ``H(lambda)`` correlations."""
    if lam >= 0.0:
        l = 0.22 + 1.57 * lam - 1.8 * lam**2
        H = 2.61 - 3.75 * lam + 5.24 * lam**2
    else:  # valid for lambda >= -0.1
        l = 0.22 + 1.402 * lam + 0.018 * lam / (lam + 0.107)
        H = 2.088 + 0.0731 / (lam + 0.14)
    return l, H


def _head_H1(H: float) -> float:
    """Entrainment shape factor ``H1`` from the shape factor ``H`` (Head)."""
    if H <= 1.6:
        return 3.3 + 0.8234 * (H - 1.1) ** -1.287
    return 3.3 + 1.5501 * (H - 0.6778) ** -3.064


def _head_H_from_H1(H1: float) -> float:
    """Invert :func:`_head_H1` to get ``H`` from ``H1``."""
    if H1 < 5.3:
        return 0.6778 + 1.1536 * (H1 - 3.3) ** -0.326
    return 1.1 + 0.86 * (H1 - 3.3) ** -0.777


def _cf_turbulent(H: float, re_theta: float) -> float:
    """Ludwieg–Tillmann turbulent skin friction (referenced to Ue)."""
    return 0.246 * 10.0 ** (-0.678 * H) * max(re_theta, 1.0) ** -0.268


# ------------------------------------------------------------------ marching


def _march_surface(name, s, x, ue, nu) -> SurfaceBL:
    """March one surface. ``s[0] = 0`` and ``ue[0] = 0`` at the stagnation point."""
    m = len(s)
    theta = np.zeros(m)
    delta_star = np.zeros(m)
    H = np.full(m, 2.59)         # Blasius value as a placeholder
    cf = np.zeros(m)

    # Thwaites needs the integral of Ue**5 from the stagnation point.
    ue5 = ue**5
    integ = np.concatenate([[0.0], np.cumsum(0.5 * (ue5[1:] + ue5[:-1]) * np.diff(s))])
    due_ds = np.gradient(ue, s)

    i_tr = -1
    i_sep = -1

    # --- laminar (Thwaites) until transition or the trailing edge ---
    i = 1
    while i < m:
        u = ue[i]
        if u <= 1e-9:
            theta[i] = theta[i - 1]
            H[i] = H[i - 1]
            i += 1
            continue
        th = np.sqrt(max(0.45 * nu / u**6 * integ[i], 0.0))
        lam = float(np.clip(th**2 * due_ds[i] / nu, -0.1, 0.25))
        l, Hi = _thwaites_l_H(lam)
        theta[i] = th
        H[i] = Hi
        delta_star[i] = Hi * th
        cf[i] = 2.0 * nu * l / (u * th)

        re_theta = u * th / nu
        re_x = u * s[i] / nu
        michel = re_x > 0 and re_theta >= 1.174 * (1.0 + 22400.0 / re_x) * re_x**0.46
        if michel or lam <= _THWAITES_SEP + 1e-9:
            i_tr = i
            break
        i += 1

    # --- turbulent (Head's entrainment method) from transition to the TE ---
    if i_tr != -1:
        # Carry theta across transition; reset the shape factor to a typical
        # attached-turbulent value.
        th = theta[i_tr]
        Ht = 1.4
        H1 = _head_H1(Ht)
        H[i_tr] = Ht
        delta_star[i_tr] = Ht * th
        cf[i_tr] = _cf_turbulent(Ht, ue[i_tr] * th / nu)

        for i in range(i_tr + 1, m):
            th, H1, Ht, sep = _head_step(
                s[i - 1], s[i], ue[i - 1], ue[i], th, H1, nu
            )
            if sep and i_sep == -1:
                i_sep = i
            theta[i] = th
            H[i] = Ht
            delta_star[i] = Ht * th
            cf[i] = _cf_turbulent(Ht, max(ue[i] * th / nu, 1.0))

    # Squire–Young: extrapolate the TE momentum thickness through the wake.
    th_te, ue_te, H_te = theta[-1], ue[-1], H[-1]
    cd = 2.0 * th_te * ue_te ** ((H_te + 5.0) / 2.0)

    x_tr = float(x[i_tr]) if i_tr != -1 else np.nan
    x_sep = float(x[i_sep]) if i_sep != -1 else np.nan
    # Separation right at the cusped TE is a method artifact; only flag it as a
    # genuine separated region when it begins forward of _SEP_X_LIMIT.
    separated = i_sep != -1 and x_sep < _SEP_X_LIMIT
    return SurfaceBL(
        name=name, s=s, x=x, ue=ue, theta=theta, delta_star=delta_star, H=H,
        cf=cf, i_transition=i_tr, x_transition=x_tr, x_separation=x_sep,
        separated=separated, cd=cd,
    )


def _head_step(s0, s1, ue0, ue1, theta, H1, nu, nsub=4):
    """Explicit sub-stepped integration of Head's two ODEs over one interval."""
    ds_total = s1 - s0
    if ds_total <= 0:
        H = _head_H_from_H1(H1)
        return theta, H1, H, H >= _TURB_SEP_H
    ds = ds_total / nsub
    due = (ue1 - ue0) / ds_total
    th = theta
    separated = False
    for k in range(nsub):
        frac = k / nsub
        ue = ue0 + (ue1 - ue0) * frac
        if ue <= 1e-9:
            continue
        H = _head_H_from_H1(H1)
        cf = _cf_turbulent(H, ue * th / nu)
        dth = cf / 2.0 - (H + 2.0) * th / ue * due
        CE = 0.0306 * max(H1 - 3.0, 1e-3) ** -0.6169
        dH1 = CE / th - H1 * due / ue - H1 * dth / th
        th = max(th + dth * ds, 1e-7)
        H1 = H1 + dH1 * ds
        if H1 <= _H1_SEP:        # turbulent separation: clamp and flag
            H1 = _H1_SEP
            separated = True
    H = _head_H_from_H1(H1)
    return th, H1, H, separated


# -------------------------------------------------------------- entry point


def _split_surfaces(sol: "Solution"):
    """Split the surface into upper/lower runs, each starting at stagnation.

    Returns ``(s, x, ue)`` arrays for the upper and lower surfaces, with the
    stagnation point prepended (``s = 0``, ``ue = 0``).
    """
    geom = sol.geom
    vt = sol.vt / sol.vinf                    # signed tangential velocity, Ue/Vinf

    # Arc length at each control point (panel centre) along the node ordering.
    s_nodes = np.concatenate([[0.0], np.cumsum(geom.length)])
    s_cp = 0.5 * (s_nodes[:-1] + s_nodes[1:])

    # The stagnation point is the sign change of vt nearest the leading edge.
    i_le = int(np.argmin(geom.xc))
    sign_changes = np.where(np.sign(vt[:-1]) != np.sign(vt[1:]))[0]
    if len(sign_changes) == 0:
        i_stag = i_le
    else:
        i_stag = int(sign_changes[np.argmin(np.abs(sign_changes - i_le))])

    # Interpolate the stagnation arc length between the two straddling points.
    v0, v1 = abs(vt[i_stag]), abs(vt[i_stag + 1])
    frac = v0 / (v0 + v1) if (v0 + v1) > 0 else 0.5
    s_stag = s_cp[i_stag] + frac * (s_cp[i_stag + 1] - s_cp[i_stag])
    x_stag = geom.xc[i_stag] + frac * (geom.xc[i_stag + 1] - geom.xc[i_stag])

    # Upper run: stagnation back toward the upper TE (decreasing index).
    up = np.arange(i_stag, -1, -1)
    s_up = np.concatenate([[0.0], s_stag - s_cp[up]])
    x_up = np.concatenate([[x_stag], geom.xc[up]])
    ue_up = np.concatenate([[0.0], np.abs(vt[up])])

    # Lower run: stagnation forward toward the lower TE (increasing index).
    lo = np.arange(i_stag + 1, geom.n)
    s_lo = np.concatenate([[0.0], s_cp[lo] - s_stag])
    x_lo = np.concatenate([[x_stag], geom.xc[lo]])
    ue_lo = np.concatenate([[0.0], np.abs(vt[lo])])

    # ``up``/``lo`` are the control-point indices for entries ``[1:]`` of each
    # run (entry 0 is the prepended stagnation point), so callers can scatter a
    # per-run quantity back onto the control points.
    return (s_up, x_up, ue_up, up), (s_lo, x_lo, ue_lo, lo)


def _friction_drag(upper: SurfaceBL, lower: SurfaceBL) -> float:
    """Integrate local skin friction into a free-stream-direction drag.

    A diagnostic subset of the total profile drag (it omits the pressure/form
    part and the wake): ``Cdf = sum cf (Ue/Vinf)^2 dx`` over both surfaces.
    """
    total = 0.0
    for surf in (upper, lower):
        # x along the run already monotonic from stagnation toward the TE.
        dx = np.abs(np.diff(surf.x))
        cf_mid = 0.5 * (surf.cf[1:] + surf.cf[:-1])
        ue_mid = 0.5 * (surf.ue[1:] + surf.ue[:-1])
        total += float(np.sum(cf_mid * ue_mid**2 * dx))
    return total


def _smooth121(a: np.ndarray, passes: int = 2) -> np.ndarray:
    """Light endpoint-preserving 1-2-1 smoothing (kills panel-to-panel noise)."""
    a = np.asarray(a, dtype=float).copy()
    for _ in range(passes):
        if len(a) < 3:
            break
        a[1:-1] = 0.25 * a[:-2] + 0.5 * a[1:-1] + 0.25 * a[2:]
    return a


def transpiration_velocity(sol: "Solution", bl: BoundaryLayer) -> np.ndarray:
    """Wall-blowing velocity that models the boundary-layer displacement effect.

    The displacement thickness is equivalent to a small normal velocity emitted
    from the surface, ``v_n = d(Ue * delta_star) / ds`` (the gradient of the
    "mass-defect" source ``m = Ue * delta_star`` along the downstream arc
    length). Returns one value per control point, in the *outward* normal
    direction and in ``Vinf`` units, ready to drop into the panel solver's
    flow-tangency right-hand side for viscous-inviscid coupling. Both runs start
    at the stagnation point where ``m = 0``.

    The mass-defect source is lightly smoothed and the result is tapered to zero
    through the cusped trailing edge (``_TE_TAPER_X`` .. TE), where ``delta*``/
    ``Ue`` are unreliable; without that, the TE spike diverges the coupling.

    This is what turns the otherwise one-way :func:`boundary_layer` estimate into
    a two-way coupling when ``panel.solve(..., couple=True)`` is used.
    """
    n = sol.geom.n
    vn = np.zeros(n)
    (_su, _xu, _ueu, i_up), (_sl, _xl, _uel, i_lo) = _split_surfaces(sol)
    for surf, idx in ((bl.upper, i_up), (bl.lower, i_lo)):
        m = _smooth121(surf.ue * surf.delta_star)   # mass-defect source, m(0)=0
        dm_ds = np.gradient(m, surf.s)              # downstream arc-length deriv.
        v = dm_ds[1:]                               # drop the prepended stagnation
        x = surf.x[1:]                              # control-point x/c for this run
        taper = np.clip((1.0 - x) / (1.0 - _TE_TAPER_X), 0.0, 1.0)
        vn[idx] = v * taper
    return np.where(np.isfinite(vn), vn, 0.0)


def boundary_layer(sol: "Solution", re: float = 1e6) -> BoundaryLayer:
    """Estimate the viscous profile drag for a solved airfoil at Reynolds ``re``.

    Parameters
    ----------
    sol:
        An inviscid :class:`~aerosim.panel.Solution`.
    re:
        Chord-based Reynolds number ``Vinf c / nu``.

    Returns
    -------
    BoundaryLayer
        Per-surface boundary-layer state plus the total profile-drag estimate.
    """
    if re <= 0:
        raise ValueError(f"Reynolds number must be positive, got {re}")
    nu = 1.0 / re  # non-dimensional kinematic viscosity (c = Vinf = 1)

    (s_up, x_up, ue_up, _), (s_lo, x_lo, ue_lo, _) = _split_surfaces(sol)
    upper = _march_surface("upper", s_up, x_up, ue_up, nu)
    lower = _march_surface("lower", s_lo, x_lo, ue_lo, nu)

    cd = upper.cd + lower.cd
    cd_friction = _friction_drag(upper, lower)
    return BoundaryLayer(re=re, upper=upper, lower=lower, cd=cd, cd_friction=cd_friction)
