"""Grid."""

import numpy as np
from ndgrid._ndgridrs import lib as _lib

_dtypes = {
    0: np.float32,
    1: np.float64,
}
_ctypes = {
    np.float32: "float",
    np.float64: "double",
}


class Grid(object):
    """Grid."""

    def __init__(self, rs_grid):
        """Initialise."""
        self._rs_grid = rs_grid

    def __del__(self):
        """Delete."""
        _lib.grid_free_grid(self._rs_grid)

    @property
    def dtype(self):
        """Data type."""
        return _dtypes[_lib.grid_dtype(self._rs_grid)]

    @property
    def _ctype(self):
        """C data type."""
        return _ctypes[self.dtype]

    @property
    def topology_dim(self):
        """Topology dimension."""
        return _lib.grid_tdim(self._rs_grid)

    @property
    def geometry_dim(self):
        """Geometry dimension."""
        return _lib.grid_gdim(self._rs_grid)
