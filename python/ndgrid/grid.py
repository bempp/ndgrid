"""Grid."""

from ndelement.reference_cell import ReferenceCellType
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


class Entity(object):
    """Entity."""

    def __init__(self, rs_entity):
        """Initialise."""
        self._rs_entity = rs_entity

    def __del__(self):
        """Delete."""
        _lib.grid_free_entity(self._rs_entity)

    @property
    def dtype(self):
        """Data type."""
        return _dtypes[_lib.entity_dtype(self._rs_entity)]

    @property
    def _ctype(self):
        """C data type."""
        return _ctypes[self.dtype]

    @property
    def local_index(self):
        """The local index of the entity."""
        return _lib.entity_local_index(self._rs_entity)

    @property
    def global_index(self):
        """The global index of the entity."""
        return _lib.entity_global_index(self._rs_entity)


class Grid(object):
    """Grid."""

    def __init__(self, rs_grid):
        """Initialise."""
        self._rs_grid = rs_grid

    def __del__(self):
        """Delete."""
        _lib.grid_free_grid(self._rs_grid)

    def entity_count(self, etype: ReferenceCellType) -> int:
        """Get the number of entities of the given type."""
        return _lib.grid_entity_count(self._rs_grid, etype.value)

    def entity(self, dim: int, local_index: int) -> Entity:
        """Get an entity from its local index."""
        return Entity(_lib.grid_entity(self._rs_grid, dim, local_index))

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
        return Grid(
            _lib.single_element_grid_new_from_raw_data_f64(
                _ffi.cast("double*", coordinates.ctypes.data),
                coordinates.shape[0],
                coordinates.shape[1],
                _ffi.cast("uintptr_t*", cells.ctypes.data),
                cells.shape[0],
                cell_type.value,
                geometry_degree,
            )
        )
    elif dtype == np.float32:
        return Grid(
            _lib.single_element_grid_new_from_raw_data_f32(
                _ffi.cast("float*", coordinates.ctypes.data),
                coordinates.shape[0],
                coordinates.shape[1],
                _ffi.cast("uintptr_t*", cells.ctypes.data),
                cells.shape[0],
                cell_type.value,
                geometry_degree,
            )
        )
    else:
        raise TypeError(f"Unsupported dtype: {dtype}")
