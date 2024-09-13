"""Grid."""

import typing
from ndelement.reference_cell import ReferenceCellType
from ndgrid.ownership import Owned, Ghost
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
    def _ctype(self) -> str:
        """C data type."""
        return _ctypes[self.dtype]

    @property
    def local_index(self) -> int:
        """The local index of the entity."""
        return _lib.entity_local_index(self._rs_entity)

    @property
    def global_index(self) -> int:
        """The global index of the entity."""
        return _lib.entity_global_index(self._rs_entity)

    # TODO: test
    @property
    def id(self) -> typing.Optional[int]:
        """The id of the entity."""
        if _lib.entity_has_id(self._rs_entity):
            return _lib.entity_id(self._rs_entity)
        else:
            return None

    # TODO: test
    @property
    def entity_type(self) -> ReferenceCellType:
        """The type of the entity."""
        return ReferenceCellType(_lib.entity_entity_type(self._rs_entity))

    # TODO: test
    @property
    def ownership(self) -> typing.Union[Owned, Ghost]:
        """The type of the entity."""
        if _lib.entity_is_owned(self._rs_entity):
            return Owned()
        else:
            return Ghost(
                _lib.entity_ownership_process(self._rs_entity),
                _lib.entity_ownership_index(self._rs_entity),
            )


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

    # TODO: test
    def entity_from_id(self, dim: int, id: int) -> Entity:
        """Get an entity from its insertion id."""
        return Entity(_lib.grid_entity_from_id(self._rs_grid, dim, id))

    # TODO: test
    def entity_types(self, dim: int) -> typing.List[ReferenceCellType]:
        """Get the entity types of the given dimension."""
        types = np.empty(_lib.entity_types_size(self._rs_grid, dim), dtype=np.uint8)
        _lib.entity_types(self._rs_grid, dim, _ffi.cast("uint8_t* ", types.ctypes.data))
        return [ReferenceCellType(i) for i in types]

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
