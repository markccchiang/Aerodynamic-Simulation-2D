"""Interactive 2D airfoil flow explorer.

Drag the sliders to change the angle of attack and the NACA 4-digit shape and
watch the streamlines, surface pressure, and force coefficients update live.

Run with:  uv run python src/main.py
"""

from __future__ import annotations

import argparse

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.widgets import Button, Slider

from .airfoil import naca4
from .flowfield import velocity_field
from .panel import Geometry, solve

# Background grid for streamlines (resolution vs. interactivity trade-off).
GRID_NX, GRID_NY = 160, 110
XLIM = (-0.6, 1.6)
YLIM = (-0.7, 0.7)
CP_RANGE = (-3.0, 1.0)


class Explorer:
    def __init__(self, alpha=5.0, m=2, p=4, t=12, n_panels=160):
        self.n_panels = n_panels
        self.xs = np.linspace(*XLIM, GRID_NX)
        self.ys = np.linspace(*YLIM, GRID_NY)

        self.fig = plt.figure(figsize=(13, 7))
        self.fig.canvas.manager.set_window_title("aerosim - airfoil flow explorer")
        gs = self.fig.add_gridspec(
            1, 2, width_ratios=[2.1, 1.0], left=0.06, right=0.97,
            top=0.93, bottom=0.30, wspace=0.22,
        )
        self.ax_flow = self.fig.add_subplot(gs[0, 0])
        self.ax_cp = self.fig.add_subplot(gs[0, 1])

        self._make_sliders()
        self._colorbar = None
        self.update(None)

    # ------------------------------------------------------------------ UI
    def _make_sliders(self):
        axc = "#dfe6ef"
        s_alpha = self.fig.add_axes([0.10, 0.18, 0.55, 0.03], facecolor=axc)
        s_m = self.fig.add_axes([0.10, 0.13, 0.55, 0.03], facecolor=axc)
        s_p = self.fig.add_axes([0.10, 0.08, 0.55, 0.03], facecolor=axc)
        s_t = self.fig.add_axes([0.10, 0.03, 0.55, 0.03], facecolor=axc)

        self.sl_alpha = Slider(s_alpha, "Angle of attack (deg)", -15, 15,
                               valinit=5.0, valstep=0.5)
        self.sl_m = Slider(s_m, "Max camber  M (%)", 0, 9, valinit=2, valstep=1)
        self.sl_p = Slider(s_p, "Camber pos.  P (x/10)", 0, 9, valinit=4, valstep=1)
        self.sl_t = Slider(s_t, "Thickness  XX (%)", 4, 30, valinit=12, valstep=1)

        for sl in (self.sl_alpha, self.sl_m, self.sl_p, self.sl_t):
            sl.on_changed(self.update)

        ax_reset = self.fig.add_axes([0.74, 0.04, 0.1, 0.05])
        self.btn_reset = Button(ax_reset, "Reset")
        self.btn_reset.on_clicked(self._reset)

    def _reset(self, _evt):
        for sl in (self.sl_alpha, self.sl_m, self.sl_p, self.sl_t):
            sl.reset()

    # -------------------------------------------------------------- compute
    def _code(self):
        m, p, t = int(self.sl_m.val), int(self.sl_p.val), int(self.sl_t.val)
        if m > 0 and p == 0:   # camber needs a valid (non-zero) position
            p = 1
            self.sl_p.eventson = False
            self.sl_p.set_val(1)
            self.sl_p.eventson = True
        return f"{m}{p}{t:02d}"

    # --------------------------------------------------------------- render
    def update(self, _val):
        code = self._code()
        alpha = float(self.sl_alpha.val)

        geom = Geometry(*naca4(code, n_panels=self.n_panels))
        sol = solve(geom, alpha)
        X, Y, U, V = velocity_field(sol, self.xs, self.ys)
        speed = np.hypot(U, V)
        cp_field = 1.0 - speed**2

        # ---- flow field ----
        self.ax_flow.clear()
        self.ax_flow.set_facecolor("#0b1021")
        strm = self.ax_flow.streamplot(
            X, Y, U, V, color=cp_field, cmap="turbo",
            density=1.4, linewidth=0.8, arrowsize=0.7,
            norm=plt.Normalize(*CP_RANGE),
        )
        self.ax_flow.fill(geom.x, geom.y, color="#e8eaed", zorder=5,
                          edgecolor="#222", linewidth=1.2)

        if self._colorbar is None:
            self._colorbar = self.fig.colorbar(
                strm.lines, ax=self.ax_flow, fraction=0.046, pad=0.02
            )
            self._colorbar.set_label("Pressure coefficient $C_p$")

        self.ax_flow.set_xlim(*XLIM)
        self.ax_flow.set_ylim(*YLIM)
        self.ax_flow.set_aspect("equal")
        self.ax_flow.set_title(
            f"NACA {code}    α = {alpha:.1f}°    "
            f"$C_l$ = {sol.cl:+.3f}    $C_m$ = {sol.cm_qc:+.3f}",
            fontsize=12,
        )
        self.ax_flow.set_xlabel("x / c")
        self.ax_flow.set_ylabel("y / c")

        # ---- surface pressure distribution ----
        self.ax_cp.clear()
        half = geom.n // 2
        xc = geom.xc
        # Node ordering is upper(TE->LE) then lower(LE->TE).
        self.ax_cp.plot(xc[:half], sol.cp[:half], color="#c0392b", lw=1.6,
                        label="upper")
        self.ax_cp.plot(xc[half:], sol.cp[half:], color="#2471a3", lw=1.6,
                        label="lower")
        self.ax_cp.axhline(0, color="0.7", lw=0.8)
        self.ax_cp.invert_yaxis()          # suction (negative Cp) plotted up
        self.ax_cp.set_xlim(-0.02, 1.02)
        self.ax_cp.set_ylim(CP_RANGE[1] + 0.5, CP_RANGE[0] - 1.0)
        self.ax_cp.set_xlabel("x / c")
        self.ax_cp.set_ylabel("$C_p$")
        self.ax_cp.set_title("Surface pressure", fontsize=11)
        self.ax_cp.legend(loc="lower right", fontsize=9)
        self.ax_cp.grid(alpha=0.25)

        self.fig.canvas.draw_idle()


def main():
    ap = argparse.ArgumentParser(description="Interactive airfoil flow explorer")
    ap.add_argument("--panels", type=int, default=160, help="number of panels")
    args = ap.parse_args()
    Explorer(n_panels=args.panels)
    plt.show()


if __name__ == "__main__":
    main()
