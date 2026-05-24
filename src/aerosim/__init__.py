"""aerosim: a 2D airfoil potential-flow (panel method) simulator."""

from .airfoil import naca4
from .panel import Geometry, solve, Solution
from .flowfield import velocity_field

__all__ = ["naca4", "Geometry", "solve", "Solution", "velocity_field"]
