# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

`aerosim` is a 2D airfoil **potential-flow** simulator: a Hess–Smith panel-method
solver plus an interactive matplotlib UI. It is **inviscid** by construction — it
predicts lift and surface pressure well at small angles but has **no drag and no
stall**. A near-zero `Cd` is therefore a *correctness check* (d'Alembert's
paradox), not a bug to fix.

## Commands

Dependencies are managed with `uv` (use it, not pip/venv directly).

```bash
uv run python src/main.py        # launch the interactive explorer (needs a display)
uv run python src/validate.py    # run the validation suite (the de-facto test suite)
uv add <package>                 # add a dependency
```

Headless UI smoke test (no display, e.g. CI or this agent's sandbox):

```bash
PYTHONPATH=src MPLBACKEND=Agg uv run python -c "from aerosim.app import Explorer; Explorer().fig.savefig('snapshot.png')"
```

There is no pytest/lint setup. `src/validate.py` is the test suite: a script of
`check()` assertions against analytical/reference aerodynamics. **Run it after any
change to `panel.py` or `airfoil.py`** — it is the safety net for the physics, and
it runs in well under a second.

## Import / layout note (easy to trip on)

Code lives under `src/` with **no `[build-system]`** in `pyproject.toml`, so the
package is *not* installed. Imports like `from aerosim ...` only resolve because
running a script puts its own directory (`src/`) on `sys.path`. Consequences:

- Run scripts **by path** (`uv run python src/main.py`), not as modules.
- `uv run python -m aerosim.app` from the repo root will **fail**. Use
  `PYTHONPATH=src` if you need module/`-m` execution or to import from elsewhere.

## Architecture

The solve pipeline is a straight line; the data contract is the `Solution`
dataclass:

```
naca4(code)  ->  Geometry  ->  solve(geom, alpha)  ->  Solution
                                                          |
                                          velocity_field(sol, xs, ys) -> grid for streamlines
```

- `src/aerosim/airfoil.py` — `naca4()` builds NACA 4-digit nodes with cosine
  spacing. **Node ordering is TE → upper surface → LE → lower surface → TE.**
  Several conventions below depend on this ordering.
- `src/aerosim/panel.py` — the physics core. `Geometry` precomputes panel
  geometry; `induced()` is the shared influence-coefficient engine; `solve()`
  assembles and solves the linear system.
- `src/aerosim/flowfield.py` — reconstructs the velocity field on a grid from a
  `Solution` (reuses `induced()`), masking points inside the body with NaN.
- `src/aerosim/app.py` — the *only* matplotlib-UI module. Everything else is
  UI-agnostic, so the solver can be driven from a script, notebook, or future web
  front-end. The `Explorer` rebuilds geometry, re-solves, and redraws on every
  slider change.

### The method (panel.py)

Each panel carries a constant-strength **source** (one unknown per panel) plus a
**single vortex strength `gamma` shared by all panels**. That is `N+1` unknowns,
closed by `N` flow-tangency conditions + **1 Kutta condition**. The Kutta row
enforces equal-and-opposite tangential velocity on the **first and last panels**
(the two trailing-edge panels — hence the node ordering matters).

`induced()` returns four `(M, N)` influence matrices used by **both** the solver
(field points = control points) and `flowfield` (field points = grid). It is the
single source of truth for the singularity math — change formulas here only.

### Sign/orientation conventions (these are subtle and were debugged carefully)

- **`_exterior_side`** is auto-detected by a point-in-polygon test, because the
  diagonal (self-influence) terms in `solve()` depend on which side of each panel
  the flow is on. If you change the node ordering or geometry source, re-verify
  this and the diagonal terms.
- **Lift** is taken from the surface-pressure integral (`cl`). The
  Kutta–Joukowski cross-check in `validate.py` uses `cl = -2·gamma·perimeter`
  (the minus sign is the CCW-positive circulation convention).
- **Moment** `cm_qc` uses the aerodynamic convention (**positive = nose-up**),
  which is the negative of the raw CCW math moment.
- **Lift-curve slope** is `~2π(1 + 0.77·t/c)` per radian, not plain `2π` — the
  `2π` flat-plate value is only the thin/zero-thickness limit. `validate.py`
  encodes the thickness-corrected expectation.

### Gotchas

- NumPy 2.x: use `np.ptp(arr)`, not `arr.ptp()`.
