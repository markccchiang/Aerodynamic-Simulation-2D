"""Load airfoil geometry from ``.dat`` coordinate files.

The solver wants nodes ordered **TE -> upper -> LE -> lower -> TE** and scaled to
unit chord (the same convention :func:`aerosim.airfoil.naca4` produces). This
module parses the two common ``.dat`` layouts, normalises to that convention,
and re-panels the points onto a cosine distribution so the panel method gets
consistent resolution regardless of how the source file was sampled.

Supported formats (auto-detected):

* **Selig** — an optional title line, then one ``x y`` pair per line running
  TE -> upper -> LE -> lower -> TE. Already in our order.
* **Lednicer** — a title line, then a count line ``n_upper n_lower``, then the
  upper surface (LE -> TE) and the lower surface (LE -> TE) as two blocks.

Public helpers: :func:`parse_dat` (text -> raw coords + name), :func:`repanel`
(cosine re-sampling), :func:`airfoil_from_dat` (path/text -> ready coords),
plus :func:`format_selig` / :func:`format_lednicer` writers and the bundled
:func:`list_samples` / :func:`load_sample_text`.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline

SAMPLE_DIR = Path(__file__).parent / "airfoils"


# ----------------------------------------------------------------- parsing


def _numeric(line: str):
    """Return the first two floats on a line, or None if it isn't a data row."""
    parts = line.replace(",", " ").split()
    if len(parts) < 2:
        return None
    try:
        return float(parts[0]), float(parts[1])
    except ValueError:
        return None


def parse_dat(text: str):
    """Parse ``.dat`` text into ``(x, y, name)`` — unit chord, TE->u->LE->l->TE.

    Auto-detects Selig vs Lednicer. ``name`` is the title line, or ``""``.
    """
    lines = text.splitlines()
    name = ""
    rows = []  # (value_pair, was_first_nonblank)
    for ln in lines:
        if not ln.strip():
            continue
        pair = _numeric(ln)
        if pair is None:
            if not name:  # first non-numeric, non-blank line is the title
                name = ln.strip()
            continue
        rows.append(pair)

    if len(rows) < 5:
        raise ValueError("could not find enough coordinate points in .dat data")

    first = rows[0]
    is_lednicer = (
        first[0] > 1.5 and first[1] > 1.5
        and abs(first[0] - round(first[0])) < 1e-3
        and abs(first[1] - round(first[1])) < 1e-3
    )

    if is_lednicer:
        n_up, n_low = int(round(first[0])), int(round(first[1]))
        body = rows[1:]
        if len(body) < n_up + n_low:
            raise ValueError("Lednicer header counts exceed the points provided")
        upper = body[:n_up]               # LE -> TE
        lower = body[n_up : n_up + n_low]  # LE -> TE
        pts = upper[::-1] + lower[1:]      # TE -> LE  then  LE -> TE (drop dup LE)
    else:
        pts = rows                         # Selig: already TE -> u -> LE -> l -> TE

    x = np.array([p[0] for p in pts], dtype=float)
    y = np.array([p[1] for p in pts], dtype=float)

    # Drop consecutive duplicate nodes (avoids zero-length panels).
    keep = np.concatenate([[True], (np.diff(x) != 0) | (np.diff(y) != 0)])
    x, y = x[keep], y[keep]

    # Ensure the upper surface comes first (some files run TE->lower->LE->upper).
    half = len(x) // 2
    if y[:half].mean() < y[half:].mean():
        x, y = x[::-1], y[::-1]

    return _normalize(x, y) + (name,)


def _normalize(x, y):
    """Translate the leading edge to x=0 and scale to unit chord (uniform)."""
    chord = float(np.ptp(x))
    if chord <= 0:
        raise ValueError("degenerate airfoil: zero chord")
    x0 = float(x.min())
    return (x - x0) / chord, y / chord


# --------------------------------------------------------------- re-paneling


def repanel(x, y, n_panels: int = 160):
    """Re-sample coordinates onto ``n_panels`` cosine-spaced panels.

    Splits at the leading edge (minimum x), fits an arc-length cubic spline to
    each surface, and resamples with cosine clustering at the leading and
    trailing edges. Falls back to the input points if a surface is too sparse
    to spline reliably.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    le = int(np.argmin(x))
    if le < 4 or (len(x) - 1 - le) < 4:
        return x, y  # too few points on a surface — keep as-is

    n_half = n_panels // 2
    t_cos = 0.5 * (1.0 - np.cos(np.linspace(0.0, np.pi, n_half + 1)))  # clustered ends

    def resample(xs, ys):
        d = np.hypot(np.diff(xs), np.diff(ys))
        s = np.concatenate([[0.0], np.cumsum(d)])
        if s[-1] == 0:
            return xs, ys
        t = s / s[-1]
        keep = np.concatenate([[True], np.diff(t) > 1e-12])  # strictly increasing
        csx = CubicSpline(t[keep], xs[keep])
        csy = CubicSpline(t[keep], ys[keep])
        return csx(t_cos), csy(t_cos)

    xu, yu = resample(x[: le + 1], y[: le + 1])  # TE -> LE
    xl, yl = resample(x[le:], y[le:])            # LE -> TE
    return np.concatenate([xu, xl[1:]]), np.concatenate([yu, yl[1:]])


def airfoil_from_dat(source, n_panels: int = 160):
    """Path or ``.dat`` text -> cosine re-paneled ``(x, y)`` ready for Geometry."""
    text = _as_text(source)
    x, y, _ = parse_dat(text)
    return repanel(x, y, n_panels)


def _as_text(source) -> str:
    if isinstance(source, Path):
        return source.read_text()
    if isinstance(source, str) and "\n" not in source and Path(source).exists():
        return Path(source).read_text()
    return source


# ----------------------------------------------------------------- writers


def format_selig(x, y, name: str = "airfoil") -> str:
    """Write coordinates as a Selig ``.dat`` string (TE -> u -> LE -> l -> TE)."""
    lines = [name]
    lines += [f"{xi:.6f} {yi:.6f}" for xi, yi in zip(x, y)]
    return "\n".join(lines) + "\n"


def format_lednicer(x, y, name: str = "airfoil") -> str:
    """Write coordinates as a Lednicer ``.dat`` string (upper/lower blocks)."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    le = int(np.argmin(x))
    upper = list(zip(x[: le + 1][::-1], y[: le + 1][::-1]))  # LE -> TE
    lower = list(zip(x[le:], y[le:]))                        # LE -> TE
    lines = [name, f"   {len(upper)}.   {len(lower)}.", ""]
    lines += [f" {xi:.6f}  {yi:.6f}" for xi, yi in upper]
    lines.append("")
    lines += [f" {xi:.6f}  {yi:.6f}" for xi, yi in lower]
    return "\n".join(lines) + "\n"


# ----------------------------------------------------------------- samples


def list_samples():
    """Return ``[(key, name), ...]`` for the bundled sample airfoils."""
    out = []
    for path in sorted(SAMPLE_DIR.glob("*.dat")):
        try:
            name = next(
                (ln.strip() for ln in path.read_text().splitlines()
                 if ln.strip() and _numeric(ln) is None),
                path.stem,
            )
        except OSError:
            name = path.stem
        out.append((path.stem, name))
    return out


def load_sample_text(key: str) -> str:
    """Return the raw ``.dat`` text of a bundled sample by its key (file stem)."""
    path = SAMPLE_DIR / f"{Path(key).stem}.dat"
    if not path.exists():
        raise FileNotFoundError(f"no bundled sample airfoil named {key!r}")
    return path.read_text()
