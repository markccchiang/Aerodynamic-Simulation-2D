# How it works

← [Back to the README](../README.md)

The airfoil surface is split into ~160 flat **panels**. Each panel carries a
constant-strength **source** (one unknown per panel) plus one shared
**vortex** strength. The unknowns are found from:

- **flow tangency** — no flow passes through the surface (one equation per panel), and
- the **Kutta condition** — flow leaves the trailing edge smoothly (one equation),

giving an `(N+1) x (N+1)` linear system. Its matrix depends only on the airfoil
*geometry* — the angle of attack enters through the right-hand side alone — so it
is factorized once per airfoil and reused, which is what makes an angle-of-attack
sweep and the viscous coupling iteration cheap. From the surface velocities we get
the pressure coefficient `Cp = 1 - (V/V∞)²`, then integrate pressure for the lift
coefficient `Cl` and moment `Cm`.

This is *inviscid potential flow*: it predicts lift and pressure distributions
well at small angles, but the panel solve itself has no viscosity, so its pressure
drag is **~0** (d'Alembert's paradox) and it has no stall.

## Viscous drag correction

On top of the inviscid solve, an **uncoupled boundary-layer correction** estimates
the *profile drag* (`Cd`) that potential flow misses. From the inviscid surface
velocities it marches an integral boundary-layer method along each surface —
**Thwaites** (laminar) → **Michel** transition → **Head's entrainment method**
(turbulent) — and applies the **Squire–Young** formula at the trailing edge to get
drag. It also locates transition and flags trailing-edge separation.

By default this is *one-way* (uncoupled): lift, pressure, and the d'Alembert check
are unchanged. Pass a chord-based Reynolds number to switch it on:

> [!NOTE]
> The code examples assume `src/` is on the import path — see
> [Using the solver from Python](using-python.md).

```python
from aerosim import naca4, Geometry, solve
sol = solve(Geometry(*naca4("0012")), alpha_deg=4.0, re=1e6)
print(sol.cd_visc)        # ~0.009 profile drag; sol.cd stays ~0 (d'Alembert)
```

In the web UI, the **Reynolds** slider drives it live: `Cd` shows in the
Coefficients panel and transition points are marked on the pressure plot.

## Two-way viscous–inviscid coupling

Add `couple=True` to feed the boundary layer **back into** the solve. The
displacement thickness is modelled as a small wall-blowing velocity
(`v_n = d(Ue·δ*)/ds`) added to the flow-tangency condition, and the solve is
iterated to convergence. Now lift and pressure **react to viscosity** — the
airfoil effectively sees a slightly thicker, decambered shape, so lift drops a
little and the suction peak softens (real "viscous decambering"):

```python
inv = solve(Geometry(*naca4("2412")), alpha_deg=5.0, re=1e6)
vis = solve(Geometry(*naca4("2412")), alpha_deg=5.0, re=1e6, couple=True)
print(inv.cl, vis.cl)          # e.g. 0.857 -> 0.773
print(vis.coupled, vis.n_iter, vis.converged)   # True 13 True
```

This is **direct** coupling: robust for attached and mildly separated flow, and it
reports `Solution.converged`. It is *not* a post-stall model — through massive
separation the direct iteration won't converge (that needs a semi-inverse scheme).
In the web UI, tick **Viscous coupling** to drive it live.

For the full derivation of all three pieces, see the Sphinx theory docs;
[the README](../README.md#building-the-sphinx-docs) shows how to build them.
