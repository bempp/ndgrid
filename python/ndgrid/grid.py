"""Grid."""

from ndelement import ReferenceCellType
import numpy as np
import numpy.typing as npt
from ndgrid._ndgridrs import lib as _lib, ffi as _ffi

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


def from_raw_data(
    coordinates: npt.NDArray[np.floating],
    cells: npt.NDArray[int],
    cell_type: ReferenceCellType,
    geometry_degree: int,
) -> Grid:
    dtype = coordinates.dtype
    if dtype == np.float64:
        return Grid(_lib.single_element_grid_new_from_raw_data_f64(
            _ffi.cast("double*", coordinates.ctypes.data),
            coordinates.shape[0],
            coordinates.shape[1],
            _ffi.cast("uintptr_t*", cells.ctypes.data),
            cells.shape[0],
            cell_type.value,
            geometry_degree,
        ))
    elif dtype == np.float32:
        return Grid(_lib.single_element_grid_new_from_raw_data_f32(
            _ffi.cast("float*", coordinates.ctypes.data),
            coordinates.shape[0],
            coordinates.shape[1],
            _ffi.cast("uintptr_t*", cells.ctypes.data),
            cells.shape[0],
            cell_type.value,
            geometry_degree,
        ))
    else:
        raise TypeError(f"Unsupported dtype: {dtype}")
