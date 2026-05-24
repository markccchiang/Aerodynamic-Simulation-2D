"""NACA 4-digit airfoil geometry and paneling.

A NACA 4-digit code ``MPXX`` encodes:
    M  -> maximum camber as a percentage of chord       (code[0])
    P  -> location of maximum camber in tenths of chord  (code[1])
    XX -> maximum thickness as a percentage of chord     (code[2:])

e.g. "2412" = 2% camber at 0.4 chord, 12% thick;  "0012" = symmetric, 12% thick.
"""

from __future__ import annotations

import numpy as np


def naca4(code: str, n_panels: int = 160, closed_te: bool = True):
    """Return node coordinates for a paneled NACA 4-digit airfoil.

    Parameters
    ----------
    code:
        4-digit NACA designation, e.g. ``"2412"``.
    n_panels:
        Number of panels around the airfoil (use an even number).
    closed_te:
        If True use the thickness coefficient that closes the trailing edge
        (better conditioned for the panel method).

    Returns
    -------
    (x, y):
        Node coordinates, each of length ``n_panels + 1``. Nodes are ordered
        from the upper trailing edge, forward to the leading edge, then back
        along the lower surface to the trailing edge, forming a closed loop.
    """
    code = code.strip()
    if len(code) != 4 or not code.isdigit():
        raise ValueError(f"NACA 4-digit code must be 4 digits, got {code!r}")

    m = int(code[0]) / 100.0       # max camber
    p = int(code[1]) / 10.0        # location of max camber
    t = int(code[2:]) / 100.0      # max thickness

    n_half = n_panels // 2
    # Cosine spacing clusters nodes near the leading and trailing edges,
    # where the geometry and pressure gradients change fastest.
    beta = np.linspace(0.0, np.pi, n_half + 1)
    xc = 0.5 * (1.0 - np.cos(beta))            # 0 -> 1

    a4 = -0.1036 if closed_te else -0.1015
    yt = 5 * t * (
        0.2969 * np.sqrt(xc)
        - 0.1260 * xc
        - 0.3516 * xc**2
        + 0.2843 * xc**3
        + a4 * xc**4
    )

    if m > 0 and p > 0:
        yc = np.where(
            xc < p,
            m / p**2 * (2 * p * xc - xc**2),
            m / (1 - p) ** 2 * ((1 - 2 * p) + 2 * p * xc - xc**2),
        )
        dyc = np.where(
            xc < p,
            2 * m / p**2 * (p - xc),
            2 * m / (1 - p) ** 2 * (p - xc),
        )
    else:  # symmetric airfoil
        yc = np.zeros_like(xc)
        dyc = np.zeros_like(xc)

    phi = np.arctan(dyc)
    xu = xc - yt * np.sin(phi)
    yu = yc + yt * np.cos(phi)
    xl = xc + yt * np.sin(phi)
    yl = yc - yt * np.cos(phi)

    # Upper surface reversed (TE -> LE) then lower surface (LE -> TE).
    x = np.concatenate([xu[::-1], xl[1:]])
    y = np.concatenate([yu[::-1], yl[1:]])
    return x, y
