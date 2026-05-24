"""aerosim: a 2D airfoil potential-flow (panel method) simulator."""

from .airfoil import naca4
from .airfoil_io import parse_dat, repanel, airfoil_from_dat, list_samples, load_sample_text
from .panel import Geometry, solve, Solution
from .flowfield import velocity_field, streamlines_from_grid
from .boundary_layer import boundary_layer, BoundaryLayer, SurfaceBL

__all__ = [
    "naca4", "Geometry", "solve", "Solution",
    "parse_dat", "repanel", "airfoil_from_dat", "list_samples", "load_sample_text",
    "velocity_field", "streamlines_from_grid",
    "boundary_layer", "BoundaryLayer", "SurfaceBL",
]
