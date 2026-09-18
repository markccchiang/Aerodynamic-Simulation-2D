# Development

← [Back to the README](../README.md)

## Project layout

| File | Responsibility |
|------|----------------|
| `src/aerosim/airfoil.py`  | NACA 4-digit geometry + cosine-spaced paneling |
| `src/aerosim/airfoil_io.py` | load `.dat` files (Selig/Lednicer), normalize + re-panel |
| `src/aerosim/airfoils/`   | bundled sample airfoil `.dat` files |
| `src/aerosim/panel.py`    | Hess–Smith solver — influence coefficients, Kutta, `Cp`/`Cl`/`Cm` |
| `src/aerosim/flowfield.py`| velocity field on a grid; matplotlib-free streamline integrator |
| `src/aerosim/boundary_layer.py` | viscous BL correction — profile drag, transition, separation; transpiration for two-way coupling |
| `src/aerosim/webapp.py`   | FastAPI backend (JSON API) for the web UI |
| `src/aerosim/static/`     | single-page Plotly.js front-end (`index.html`, `app.js`, `styles.css`) |
| `src/web.py`              | entry point — web server (uvicorn) |
| `src/validate.py`         | checks against thin-airfoil theory & known results |
| `docs/theory/`            | Sphinx write-up of the theory, with solver-generated figures |
| `guide/`                  | these Markdown guides |
| `gallery/`                | screenshots and diagrams used by the README and the guides |

The physics core is UI-agnostic — you can drive it from a script or a notebook,
not just the web UI (see [Using the solver from Python](../README.md#using-the-solver-from-python)).

## Possible next steps

- **Semi-inverse coupling** to push the viscous solve through separation toward
  real post-stall behaviour — the current `couple=True` is *direct* coupling, so
  it captures viscous decambering but not massive separation.
- A drag breakdown (friction vs form vs induced) and a `Cl`/`Cd` (L/D) readout
  on the polar.
