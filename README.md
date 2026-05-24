# aerosim

An interactive **2D airfoil aerodynamics simulator** built from scratch in Python.
It solves the incompressible potential-flow field around a NACA airfoil with a
**Hess–Smith panel method** and lets you explore lift, pressure, and streamlines
in real time with sliders.

![preview](snapshot.png)

## Quick start

```bash
uv run python src/main.py        # launch the interactive desktop explorer
uv run python src/web.py         # launch the web UI -> http://127.0.0.1:8000
uv run python src/validate.py    # run the physics validation suite
```

Drag the sliders to change the angle of attack, Reynolds number, and the NACA
4-digit shape (`MPXX` = camber %, camber position, thickness %). The flow field,
surface pressure plot, and the lift/drag/moment coefficients all update live.

## Web UI

The simulator also runs in the browser. `src/web.py` starts a small **FastAPI**
backend that drives the same solver and returns JSON (geometry, a `Cp` field,
streamlines, surface pressure, and coefficients); a single **Plotly.js** page
renders it with live sliders. No build step — open the printed URL.

![web preview](web-preview.png)

The desktop matplotlib explorer and the web UI are two thin front-ends over the
same UI-agnostic solver core.

## How it works

The airfoil surface is split into ~160 flat **panels**. Each panel carries a
constant-strength **source** (one unknown per panel) plus one shared
**vortex** strength. The unknowns are found from:

- **flow tangency** — no flow passes through the surface (one equation per panel), and
- the **Kutta condition** — flow leaves the trailing edge smoothly (one equation),

giving an `(N+1) x (N+1)` linear system. From the surface velocities we get the
pressure coefficient `Cp = 1 - (V/V∞)²`, then integrate pressure for the lift
coefficient `Cl` and moment `Cm`.

This is *inviscid potential flow*: it predicts lift and pressure distributions
well at small angles, but the panel solve itself has no viscosity, so its pressure
drag is **~0** (d'Alembert's paradox) and it has no stall.

### Viscous drag correction

On top of the inviscid solve, an **uncoupled boundary-layer correction** estimates
the *profile drag* (`Cd`) that potential flow misses. From the inviscid surface
velocities it marches an integral boundary-layer method along each surface —
**Thwaites** (laminar) → **Michel** transition → **Head's entrainment method**
(turbulent) — and applies the **Squire–Young** formula at the trailing edge to get
drag. It also locates transition and flags trailing-edge separation.

It is *one-way*: lift, pressure, and the d'Alembert check are unchanged. Pass a
chord-based Reynolds number to switch it on:

```python
from aerosim import naca4, Geometry, solve
sol = solve(Geometry(*naca4("0012")), alpha_deg=4.0, re=1e6)
print(sol.cd_visc)        # ~0.009 profile drag; sol.cd stays ~0 (d'Alembert)
```

In the explorer, the **log10(Reynolds)** slider drives it live: `Cd` shows in the
title and transition points are marked on the pressure plot.

## Project layout

| File | Responsibility |
|------|----------------|
| `src/aerosim/airfoil.py`  | NACA 4-digit geometry + cosine-spaced paneling |
| `src/aerosim/panel.py`    | Hess–Smith solver — influence coefficients, Kutta, `Cp`/`Cl`/`Cm` |
| `src/aerosim/flowfield.py`| velocity field on a grid; matplotlib-free streamline integrator |
| `src/aerosim/boundary_layer.py` | uncoupled viscous BL correction — profile drag, transition, separation |
| `src/aerosim/app.py`      | interactive matplotlib desktop UI |
| `src/aerosim/webapp.py`   | FastAPI backend (JSON API) for the web UI |
| `src/aerosim/static/index.html` | single-page Plotly.js web front-end |
| `src/main.py`             | entry point — desktop explorer |
| `src/web.py`              | entry point — web server (uvicorn) |
| `src/validate.py`         | checks against thin-airfoil theory & known results |

The physics core (`panel.py`) is UI-agnostic — you can drive it from a script,
a notebook, or a future web front-end.

## Validation

This is the project's test suite — **run `uv run python src/validate.py` after any change to the physics core** (`src/aerosim/panel.py`, `airfoil.py`, or `boundary_layer.py`); it runs in well under a second.

`validate.py` confirms the solver against analytical / reference results:

- symmetric airfoil produces zero lift at zero incidence;
- lift-curve slope ≈ `2π(1 + 0.77·t/c)` per radian;
- cambered NACA 2412 gives `Cl ≈ 0.25` at α = 0 and a nose-down `Cm`;
- pressure drag ≈ 0 (d'Alembert);
- peak surface `Cp ≈ 1` at the stagnation point;
- lift from the pressure integral matches the Kutta–Joukowski circulation.

For the viscous correction it also checks:

- profile `Cd` of NACA 0012 ≈ published data at Re = 10⁶;
- `Cd` falls with Reynolds number and rises with angle of attack;
- transition advances toward the leading edge as Reynolds number rises;
- the inviscid `Cl`/`Cd` are untouched by enabling the viscous estimate;
- the flow is attached at α = 0 and separates as it approaches stall.

## Possible next steps

- **Couple** the boundary layer back into the inviscid solve (displacement-
  thickness transpiration) so lift and pressure react to viscosity and real
  stall appears — the current correction is one-way (drag only).
- Load arbitrary airfoil coordinate files (`.dat`) instead of only NACA shapes.
- Plot `Cl` vs α and the drag polar (`Cl` vs `Cd`) now that `Cd` exists.
