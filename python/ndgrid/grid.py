"""Grid."""

import typing
from ndelement.reference_cell import ReferenceCellType
from ndgrid.ownership import Owned, Ghost
import numpy as np
import numpy.typing as npt
from ndgrid._ndgridrs import lib as _lib, ffi as _ffi
from _cffi_backend import _CDataBase

_dtypes = {
    0: np.float32,
    1: np.float64,
}
_ctypes = {
    np.float32: "float",
    np.float64: "double",
}


class Topology(object):
    """Entity topology."""

    def __init__(self, rs_topology: _CDataBase, owned: bool = True):
        """Initialise."""
        self._rs_topology = rs_topology
        self._owned = owned

    def __del__(self):
        """Delete."""
        if self._owned:
            _lib.free_topology(self._rs_topology)

    def sub_entity(self, dim: int, index: int) -> int:
        """Get the index of a subentity."""
        return _lib.topology_sub_entity(self._rs_topology, dim, index)

    def sub_entities(self, dim: int) -> typing.List[int]:
        """Get points that define the geometry."""
        entities = np.empty(_lib.topology_sub_entities_size(self._rs_topology, dim), dtype=np.uintp)
        _lib.topology_sub_entities(
            self._rs_topology, dim, _ffi.cast("uintptr_t* ", entities.ctypes.data)
        )
        return [int(i) for i in entities]

    def connected_entities(self, dim: int) -> typing.List[int]:
        """Get points that define the geometry."""
        entities = np.empty(
            _lib.topology_connected_entities_size(self._rs_topology, dim), dtype=np.uintp
        )
        _lib.topology_connected_entities(
            self._rs_topology, dim, _ffi.cast("uintptr_t* ", entities.ctypes.data)
        )
        return [int(i) for i in entities]


class Geometry(object):
    """Entity geometry."""

    def __init__(self, rs_geometry: _CDataBase, dim: int, owned: bool = True):
        """Initialise."""
        self._rs_geometry = rs_geometry
        self._dim = dim
        self._owned = owned

    def __del__(self):
        """Delete."""
        if self._owned:
            _lib.free_geometry(self._rs_geometry)

    def points(self) -> npt.NDArray:
        """Get points that define the geometry."""
        pts = np.empty((self.point_count, self._dim), dtype=self.dtype)
        _lib.geometry_points(self._rs_geometry, _ffi.cast("void* ", pts.ctypes.data))
        return pts

    @property
    def dtype(self):
        """Data type."""
        return _dtypes[_lib.geometry_dtype(self._rs_geometry)]

    @property
    def _ctype(self) -> str:
        """C data type."""
        return _ctypes[self.dtype]

    @property
    def point_count(self) -> int:
        """Number of points."""
        return _lib.geometry_point_count(self._rs_geometry)

    @property
    def degree(self) -> int:
        """Degree of geometry element."""
        return _lib.geometry_degree(self._rs_geometry)


class GeometryMap(object):
    """Geometry map."""

    def __init__(self, rs_gmap: _CDataBase, owned: bool = True):
        """Initialise."""
        self._rs_gmap = rs_gmap
        self._owned = owned

    def __del__(self):
        """Delete."""
        if self._owned:
            _lib.free_geometry_map(self._rs_gmap)

    def points(self, entity_index: int, points: npt.NDArray[np.floating]):
        """Get the physical points for an entity."""
        assert points.dtype == self.dtype
        _lib.geometry_map_points(
            self._rs_gmap, entity_index, _ffi.cast("void* ", points.ctypes.data)
        )

    def jacobians(self, entity_index: int, jacobians: npt.NDArray[np.floating]):
        """Get the jacobians for an entity."""
        assert jacobians.dtype == self.dtype
        _lib.geometry_map_jacobians(
            self._rs_gmap, entity_index, _ffi.cast("void* ", jacobians.ctypes.data)
        )

    def jacobians_dets_normals(
        self,
        entity_index: int,
        jacobians: npt.NDArray[np.floating],
        jdets: npt.NDArray[np.floating],
        normals: npt.NDArray[np.floating],
    ):
        """Get the jacobians, their determinants, and the normals for an entity."""
        assert jacobians.dtype == self.dtype
        assert jdets.dtype == self.dtype
        assert normals.dtype == self.dtype
        _lib.geometry_map_jacobians_dets_normals(
            self._rs_gmap,
            entity_index,
            _ffi.cast("void* ", jacobians.ctypes.data),
            _ffi.cast("void* ", jdets.ctypes.data),
            _ffi.cast("void* ", normals.ctypes.data),
        )

    @property
    def dtype(self):
        """Data type."""
        return _dtypes[_lib.geometry_map_dtype(self._rs_gmap)]

    @property
    def _ctype(self) -> str:
        """C data type."""
        return _ctypes[self.dtype]

    @property
    def entity_topology_dimension(self) -> int:
        """The topoloical dimension of the entity being mapped."""
        return _lib.geometry_map_entity_topology_dimension(self._rs_gmap)

    @property
    def geometry_dimension(self) -> int:
        """The geometric dimension of the physical space."""
        return _lib.geometry_map_geometry_dimension(self._rs_gmap)

    @property
    def point_count(self) -> int:
        """The number of reference points that this map uses."""
        return _lib.geometry_map_point_count(self._rs_gmap)


class Entity(object):
    """Entity."""

    def __init__(self, rs_entity: _CDataBase, gdim: int, tdim: int, owned: bool = True):
        """Initialise."""
        self._rs_entity = rs_entity
        self._gdim = gdim
        self._tdim = tdim
        self._owned = owned

    def __del__(self):
        """Delete."""
        if self._owned:
            _lib.free_entity(self._rs_entity)

    @property
    def dtype(self):
        """Data type."""
        return _dtypes[_lib.entity_dtype(self._rs_entity)]

    @property
    def _ctype(self) -> str:
        """C data type."""
        return _ctypes[self.dtype]

    @property
    def topology(self) -> Topology:
        """Entity topology."""
        return Topology(_lib.entity_topology(self._rs_entity))

    @property
    def geometry(self) -> Geometry:
        """Entity geometry."""
        return Geometry(_lib.entity_geometry(self._rs_entity), self._gdim)

    @property
    def local_index(self) -> int:
        """The local index of the entity."""
        return _lib.entity_local_index(self._rs_entity)

    @property
    def global_index(self) -> int:
        """The global index of the entity."""
        return _lib.entity_global_index(self._rs_entity)

    @property
    def id(self) -> typing.Optional[int]:
        """The id of the entity."""
        if _lib.entity_has_id(self._rs_entity):
            return _lib.entity_id(self._rs_entity)
        else:
            return None

    @property
    def entity_type(self) -> ReferenceCellType:
        """The type of the entity."""
        return ReferenceCellType(_lib.entity_entity_type(self._rs_entity))

    @property
    def ownership(self) -> typing.Union[Owned, Ghost]:
        """The ownership of the entity."""
        if _lib.entity_is_owned(self._rs_entity):
            return Owned()
        else:
            return Ghost(
                _lib.entity_ownership_process(self._rs_entity),
                _lib.entity_ownership_index(self._rs_entity),
            )


class Grid(object):
    """Grid."""

    def __init__(self, rs_grid: _CDataBase, owned: bool = True):
        """Initialise."""
        self._rs_grid = rs_grid
        self._owned = owned

    def __del__(self):
        """Delete."""
        if self._owned:
            _lib.free_grid(self._rs_grid)

    def entity_count(self, etype: ReferenceCellType) -> int:
        """Get the number of entities of the given type."""
        return _lib.grid_entity_count(self._rs_grid, etype.value)

    def entity(self, dim: int, local_index: int) -> Entity:
        """Get an entity from its local index."""
        return Entity(_lib.grid_entity(self._rs_grid, dim, local_index), self.geometry_dim, dim)

    def entity_from_id(self, dim: int, id: int) -> Entity:
        """Get an entity from its insertion id."""
        return Entity(_lib.grid_entity_from_id(self._rs_grid, dim, id), self.geometry_dim, dim)

    def entity_types(self, dim: int) -> typing.List[ReferenceCellType]:
        """Get the entity types of the given dimension."""
        types = np.empty(_lib.grid_entity_types_size(self._rs_grid, dim), dtype=np.uint8)
        _lib.grid_entity_types(self._rs_grid, dim, _ffi.cast("uint8_t* ", types.ctypes.data))
        return [ReferenceCellType(i) for i in types]

    def geometry_map(
        self, entity_type: ReferenceCellType, points: npt.NDArray[np.floating]
    ) -> GeometryMap:
        """Get a geometry map for the given points on an entity."""
        assert points.dtype == self.dtype
        return GeometryMap(
            _lib.grid_geometry_map(
                self._rs_grid,
                entity_type.value,
                _ffi.cast("void* ", points.ctypes.data),
                points.shape[0],
            )
        )

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
    cells: npt.NDArray[np.int64],
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
