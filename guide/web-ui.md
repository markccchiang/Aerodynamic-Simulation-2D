# Web UI

← [Back to the README](../README.md)

`src/web.py` starts a small **FastAPI** backend that drives the solver and returns
JSON (geometry, a `Cp` field, streamlines, surface pressure, and coefficients); a
single **Plotly.js** page renders it with live sliders. No build step — open the
printed URL. Besides the NACA sliders, the **Airfoil** panel lets you pick a
bundled sample or **upload your own `.dat` file** (see
[Loading airfoils](airfoils.md)).

The web page is a thin front-end — all physics lives in the UI-agnostic solver
core, which you can also drive directly from a script or notebook (see
[Using the solver from Python](../README.md#using-the-solver-from-python)).

## Lift curve & drag polar

Below the flow view, a **Lift curve & drag polar** panel sweeps the angle of
attack (the **Compute polar** button) and plots `Cl` vs α and the drag polar
(`Cl` vs `Cd`), marking the current operating point. With **Viscous coupling**
on, it overlays the inviscid and coupled curves so the viscous decambering is
visible directly; coupled points where the iteration did not converge (it stops
converging as you approach stall) are ringed in amber rather than drawn as
settled answers.

## Animation

Between the flow view and the polar sits an **Animation** strip: choose angle of
attack, max camber M, camber position P or thickness, give it a from/to/step
range, and **Build** sweeps that one parameter and plays the result back — the
flow field, surface pressure and coefficients all follow the sweep. Frames are
solved once and cached, so scrubbing, looping and changing speed never re-solve.
The sliders are left alone while it plays (the flow plot's title and the
animation strip show which frame is on screen), and touching any of them stops
playback and returns to live mode. A loaded `.dat` airfoil can only sweep α,
since M/P/t are NACA shape digits.
