# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

`aerosim` is a 2D airfoil **potential-flow** simulator: a Hess–Smith panel-method
solver plus an interactive matplotlib UI. The panel solver is **inviscid** by
construction — it predicts lift and surface pressure well at small angles but has
**no drag and no stall**. A near-zero `Solution.cd` is therefore a *correctness
check* (d'Alembert's paradox), not a bug to fix.

On top of that core there is an **optional, uncoupled viscous boundary-layer
correction** (`boundary_layer.py`) that estimates *profile drag* the inviscid
model cannot. It reads the inviscid edge velocities and marches an integral BL
method, but does **not** feed back into the panel solve — so `cl`, `cp` and the
d'Alembert `cd` are unchanged whether or not it runs. Pass `re=` to `solve()` to
turn it on; the result lands in `Solution.bl` and `Solution.cd_visc`. Two `cd`s
therefore coexist and mean different things: `cd` (inviscid pressure drag, ≈ 0,
the d'Alembert check) vs `cd_visc` / `bl.cd` (the viscous profile-drag estimate).

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
`check()` assertions against analytical/reference aerodynamics (including a
viscous block: profile-drag magnitude vs published NACA0012 data, plus Re/α
trends, transition movement and stall onset). **Run it after any change to
`panel.py`, `airfoil.py`, or `boundary_layer.py`** — it is the safety net for the
physics, and it runs in well under a second.

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
naca4(code)  ->  Geometry  ->  solve(geom, alpha[, re])  ->  Solution
                                                              |
                                  velocity_field(sol, xs, ys) -> grid for streamlines
                                  boundary_layer(sol, re)     -> BoundaryLayer (profile drag)
```

- `src/aerosim/airfoil.py` — `naca4()` builds NACA 4-digit nodes with cosine
  spacing. **Node ordering is TE → upper surface → LE → lower surface → TE.**
  Several conventions below depend on this ordering.
- `src/aerosim/panel.py` — the physics core. `Geometry` precomputes panel
  geometry; `induced()` is the shared influence-coefficient engine; `solve()`
  assembles and solves the linear system.
- `src/aerosim/flowfield.py` — reconstructs the velocity field on a grid from a
  `Solution` (reuses `induced()`), masking points inside the body with NaN.
- `src/aerosim/boundary_layer.py` — the uncoupled viscous correction. Splits the
  surface at the stagnation point, then per surface marches **Thwaites** (laminar)
  → **Michel** transition → **Head's entrainment method** (turbulent, with
  Ludwieg–Tillmann `Cf`) and gets profile drag from the **Squire–Young** formula
  at the trailing edge. UI-agnostic, like the rest of the core. See its own
  conventions below.
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

### The viscous correction (boundary_layer.py)

- **Uncoupled by design.** It is a one-way post-process: inviscid `Solution` →
  edge velocities → BL → drag. It never modifies `sigma`/`gamma`/`cp`/`cl`. If you
  ever want stall/`Cl`-correction you'd need *viscous–inviscid coupling*
  (displacement-thickness transpiration), which is a much bigger, convergence-
  sensitive change — don't bolt it on here without re-thinking the validation.
- **Non-dimensionalisation.** Everything is in chord (`c = 1`) and `Vinf` units,
  so kinematic viscosity is simply `nu = 1/Re` and edge velocity is `Ue/Vinf`.
  Reynolds number is **chord-based**.
- **Edge velocity / stagnation.** `Ue` per surface is `|sol.vt|`; the surface is
  split at the stagnation point, found as the sign change of `vt` nearest the
  leading edge. This relies on the same TE→upper→LE→lower→TE node ordering as the
  solver.
- **The `separated` flag ignores the last ~10% of chord.** The closed/cusped TE
  makes `Ue` drop steeply over the final ~0.3% chord, which spikes the shape
  factor and trips Head's separation criterion artificially. So separation is only
  reported as genuine when it begins forward of `_SEP_X_LIMIT` (0.90 c). Don't
  "fix" the resulting near-TE separation — it's expected.
- **Validation is mostly trends, not tight tolerances.** Integral BL methods are
  approximate; `validate.py` checks one banded absolute `Cd` against published
  data and otherwise asserts physical *trends* (Cd↓ with Re, Cd↑ with α,
  transition advancing with Re, friction < profile drag). Keep new checks in that
  spirit.

### Gotchas

- NumPy 2.x: use `np.ptp(arr)`, not `arr.ptp()`.
