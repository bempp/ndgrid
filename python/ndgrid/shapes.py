"""Shapes."""

import numpy as np
import typing
from ndgrid._ndgridrs import lib as _lib
from ndgrid.grid import Grid


def regular_sphere(level: int, dtype: typing.Type[np.floating] = np.float64) -> Grid:
    """Create a regular sphere."""
    if dtype == np.float64:
        return Grid(_lib.shapes_regular_sphere_f64(level))
    elif dtype == np.float32:
        return Grid(_lib.shapes_regular_sphere_f32(level))
    else:
        raise TypeError(f"Unsupported dtype: {dtype}")
