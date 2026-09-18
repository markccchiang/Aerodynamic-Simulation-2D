# Using the solver from Python

← [Back to the README](../README.md)

## Running your own code

The `aerosim` package isn't installed, so to import it from your own script,
run from the repo root with `src/` on the import path:

```bash
PYTHONPATH=src uv run python my_script.py
```

In a notebook, `sys.path.insert(0, "src")` (from the repo root) does the same.
The code examples in the guides assume one or the other.

## A worked example

A typical `my_script.py`, from geometry to a pressure plot:

```python
from aerosim import naca4, airfoil_from_dat, Geometry, solve

# 1. Geometry: a NACA 2412 split into 160 flat panels. Nodes run
#    trailing edge -> upper surface -> leading edge -> lower surface -> trailing edge.
geom = Geometry(*naca4("2412", n_panels=160))

# 2. Inviscid solve at 5 degrees angle of attack.
sol = solve(geom, alpha_deg=5.0)
print(f"inviscid  Cl = {sol.cl:.3f}   Cm(c/4) = {sol.cm_qc:+.3f}   Cd = {sol.cd:+.5f}")

# 3. Give a chord Reynolds number and the boundary layer estimates profile drag.
visc = solve(geom, alpha_deg=5.0, re=1e6)
print(f"viscous   Cd = {visc.cd_visc:.4f}   transition x/c: "
      f"upper {visc.bl.upper.x_transition:.2f}, lower {visc.bl.lower.x_transition:.2f}")

# 4. couple=True feeds the boundary layer back into the solve.
cpl = solve(geom, alpha_deg=5.0, re=1e6, couple=True)
print(f"coupled   Cl = {cpl.cl:.3f}   converged: {cpl.converged} after {cpl.n_iter} iterations")

# 5. An angle-of-attack sweep. Reusing one Geometry means the panel matrix is
#    factorized once and every later solve is cheap.
for alpha in range(-4, 13, 4):
    s = solve(geom, alpha_deg=alpha, re=1e6)
    print(f"  alpha = {alpha:3d}   Cl = {s.cl:+.3f}   Cd = {s.cd_visc:.4f}")

# 6. Any .dat file works the same way (Selig or Lednicer layout).
x, y = airfoil_from_dat("src/aerosim/airfoils/naca4415.dat")
print(f"NACA 4415 at 5 deg: Cl = {solve(Geometry(x, y), alpha_deg=5.0).cl:.3f}")

# 7. Surface pressure is one value per panel, at the control points geom.xc, geom.yc.
import matplotlib.pyplot as plt
plt.plot(geom.xc, sol.cp)
plt.gca().invert_yaxis()                      # suction (negative Cp) plotted upward
plt.xlabel("x / c"); plt.ylabel("Cp")
plt.savefig("cp.png")
```

It prints

```text
inviscid  Cl = 0.857   Cm(c/4) = -0.061   Cd = -0.00027
viscous   Cd = 0.0108   transition x/c: upper 0.18, lower 0.90
coupled   Cl = 0.773   converged: True after 13 iterations
  alpha =  -4   Cl = -0.227   Cd = 0.0086
  alpha =   0   Cl = +0.255   Cd = 0.0077
  alpha =   4   Cl = +0.737   Cd = 0.0095
  alpha =   8   Cl = +1.215   Cd = 0.0170
  alpha =  12   Cl = +1.688   Cd = 0.0284
NACA 4415 at 5 deg: Cl = 1.141
```

and saves the surface-pressure plot to `cp.png`.

## What a solve returns

Everything a solve returns is on the `Solution` object:

| Field | What it holds |
|-------|---------------|
| `cl`, `cm_qc` | lift coefficient; moment about the quarter chord (positive nose-up) |
| `cd` | pressure drag — about zero without viscosity (d'Alembert's paradox) |
| `cp`, `vt` | surface pressure coefficient and tangential velocity, one per panel, at `geom.xc`, `geom.yc` |
| `cd_visc` | profile drag from the boundary layer (`None` unless `re=` is given) |
| `bl` | the boundary layer itself: `bl.upper` and `bl.lower` carry `x_transition`, `x_separation`, `theta`, `H`, `cf`, … |
| `coupled`, `n_iter`, `converged` | coupling status, when `couple=True` |

## Your own coordinates, and the flow field

For your own coordinates, pass a closed loop of nodes — unit chord, leading edge
at x = 0, in the order above — straight to `Geometry(x, y)`, or save them as a
`.dat` file and let `airfoil_from_dat()` normalise and re-panel them. To get the
flow off the surface, `velocity_field(sol, xs, ys)` evaluates the velocity on a
grid.
