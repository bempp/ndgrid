//! Binding for C
#![allow(missing_docs)]
#![allow(clippy::missing_safety_doc)]

use c_api_tools::cfuncs;

/// Wrapper for grid types.
#[cfuncs(name = "grid_t", create, free, unwrap)]
pub struct GridT;

/// Wrapper for entities.
#[cfuncs(name = "entity_t", create, free, unwrap)]
pub struct EntityT;

/// Wrapper for entity topology.
#[cfuncs(name = "topology_t", create, free, unwrap)]
pub struct TopologyT;

/// Wrapper for entity geometry.
#[cfuncs(name = "geometry_t", create, free, unwrap)]
pub struct GeometryT;

/// Wrapper for geometry maps.
#[cfuncs(name = "geometry_map_t", create, free, unwrap)]
pub struct GeometryMapT;

pub mod grid {
    use super::{
        entity_t_create, entity_t_unwrap, geometry_map_t_create, geometry_map_t_unwrap,
        grid_t_create, grid_t_unwrap, EntityT, GeometryMapT, GridT,
    };
    use crate::{traits::Grid, types::RealScalar, SingleElementGrid};
    use c_api_tools::{concretise_types, DType, DTypeIdentifier};
    use ndelement::{
        ciarlet::CiarletElement,
        ciarlet::LagrangeElementFamily,
        reference_cell,
        traits::{ElementFamily, FiniteElement},
        types::{Continuity, ReferenceCellType},
    };
    use std::ffi::c_void;
    use std::slice::from_raw_parts;

    #[no_mangle]
    #[allow(clippy::too_many_arguments)]
    pub extern "C" fn single_element_grid_new(
        coordinates: *const c_void,
        npts: usize,
        gdim: usize,
        cells: *const usize,
        ncells: usize,
        cell_type: ReferenceCellType,
        geometry_degree: usize,
        dtype: DType,
    ) -> *mut GridT {
        pub fn single_elemeng_grid_impl<T: RealScalar>(
            coordinates: *const T,
            npts: usize,
            gdim: usize,
            cells: *const usize,
            ncells: usize,
            cell_type: ReferenceCellType,
            geometry_degree: usize,
        ) -> *mut GridT {
            let wrapper = grid_t_create();
            let inner = unsafe { grid_t_unwrap(wrapper) }.unwrap();
            let coordinates = unsafe { from_raw_parts(coordinates, npts * gdim) };
            let points_per_cell =
                LagrangeElementFamily::<T>::new(geometry_degree, Continuity::Standard)
                    .element(cell_type)
                    .dim();
            let cells = unsafe { from_raw_parts(cells, ncells * points_per_cell) };

            *inner = Box::new(
                SingleElementGrid::<T, CiarletElement<T>>::new_from_raw_data(
                    coordinates,
                    gdim,
                    cells,
                    cell_type,
                    geometry_degree,
                ),
            );
            wrapper
        }

        match dtype {
            DType::F32 => single_elemeng_grid_impl(
                coordinates as *const f32,
                npts,
                gdim,
                cells,
                ncells,
                cell_type,
                geometry_degree,
            ),
            DType::F64 => single_elemeng_grid_impl(
                coordinates as *const f64,
                npts,
                gdim,
                cells,
                ncells,
                cell_type,
                geometry_degree,
            ),
            _ => panic!("Unsupported dtype"),
        }
    }

    #[concretise_types(
        gen_type(name = "dtype", replace_with = ["f32", "f64"]),
        field(arg = 0, name = "wrap", wrapper = "GridT", replace_with = ["SingleElementGrid<{{dtype}}, CiarletElement<{{dtype}}>>"]),
    )]
    pub fn grid_tdim<G: Grid>(grid: &G) -> usize {
        grid.topology_dim()
    }

    #[concretise_types(
        gen_type(name = "dtype", replace_with = ["f32", "f64"]),
        field(arg = 0, name = "wrap", wrapper = "GridT", replace_with = ["SingleElementGrid<{{dtype}}, CiarletElement<{{dtype}}>>"]),
    )]
    pub fn grid_gdim<G: Grid>(grid: &G) -> usize {
        grid.geometry_dim()
    }

    #[concretise_types(
        gen_type(name = "dtype", replace_with = ["f32", "f64"]),
        field(arg = 0, name = "wrap", wrapper = "GridT", replace_with = ["SingleElementGrid<{{dtype}}, CiarletElement<{{dtype}}>>"]),
    )]
    pub fn grid_entity_count<G: Grid<EntityDescriptor = ReferenceCellType>>(
        grid: &G,
        entity_type: ReferenceCellType,
    ) -> usize {
        grid.entity_count(entity_type)
    }

    #[concretise_types(
        gen_type(name = "dtype", replace_with = ["f32", "f64"]),
        field(arg = 0, name = "wrap", wrapper = "GridT", replace_with = ["SingleElementGrid<{{dtype}}, CiarletElement<{{dtype}}>>"]),
    )]
    pub fn grid_entity<G: Grid>(
        grid: &'static G,
        dim: usize,
        local_index: usize,
    ) -> *const EntityT {
        let wrapper = entity_t_create();
        let inner = unsafe { entity_t_unwrap(wrapper) }.unwrap();
        let entity = grid.entity(dim, local_index);

        if let Some(entity) = entity {
            *inner = Box::new(entity);
        }

        wrapper
    }

    #[concretise_types(
        gen_type(name = "dtype", replace_with = ["f32", "f64"]),
        field(arg = 0, name = "wrap", wrapper = "GridT", replace_with = ["SingleElementGrid<{{dtype}}, CiarletElement<{{dtype}}>>"]),
    )]
    pub fn grid_entity_from_id<G: Grid>(grid: &'static G, dim: usize, id: usize) -> *const EntityT {
        let wrapper = entity_t_create();
        let inner = unsafe { entity_t_unwrap(wrapper) }.unwrap();
        let entity = grid.entity_from_id(dim, id);

        if let Some(entity) = entity {
            *inner = Box::new(entity);
        }

        wrapper
    }

    #[concretise_types(
        gen_type(name = "dtype", replace_with = ["f32", "f64"]),
        field(arg = 0, name = "wrap", wrapper = "GridT", replace_with = ["SingleElementGrid<{{dtype}}, CiarletElement<{{dtype}}>>"]),
    )]
    pub fn grid_geometry_map<G: Grid<EntityDescriptor = ReferenceCellType>>(
        grid: &'static G,
        entity_type: ReferenceCellType,
        points: *const c_void,
        npoints: usize,
    ) -> *const GeometryMapT {
        let points = points as *const G::T;
        let wrapper = geometry_map_t_create();
        let inner = unsafe { geometry_map_t_unwrap(wrapper) }.unwrap();

        *inner = Box::new(grid.geometry_map(entity_type, unsafe {
            from_raw_parts(points, npoints * reference_cell::dim(entity_type))
        }));

        wrapper
    }

    #[concretise_types(
        gen_type(name = "dtype", replace_with = ["f32", "f64"]),
        field(arg = 0, name = "wrap", wrapper = "GridT", replace_with = ["SingleElementGrid<{{dtype}}, CiarletElement<{{dtype}}>>"]),
    )]
    pub fn grid_entity_types_size<G: Grid>(grid: &G, dim: usize) -> usize {
        grid.entity_types(dim).len()
    }

    #[concretise_types(
        gen_type(name = "dtype", replace_with = ["f32", "f64"]),
        field(arg = 0, name = "wrap", wrapper = "GridT", replace_with = ["SingleElementGrid<{{dtype}}, CiarletElement<{{dtype}}>>"]),
    )]
    pub fn grid_entity_types<G: Grid<EntityDescriptor = ReferenceCellType>>(
        grid: &G,
        dim: usize,
        types: *mut ReferenceCellType,
    ) {
        for (i, t) in grid.entity_types(dim).iter().enumerate() {
            unsafe {
                *types.add(i) = *t;
            }
        }
    }

    #[concretise_types(
        gen_type(name = "dtype", replace_with = ["f32", "f64"]),
        field(arg = 0, name = "wrap", wrapper = "GridT", replace_with = ["SingleElementGrid<{{dtype}}, CiarletElement<{{dtype}}>>"]),
    )]
    pub fn grid_dtype<T: RealScalar + DTypeIdentifier, G: Grid<T = T>>(_grid: &G) -> DType {
        <T as DTypeIdentifier>::dtype()
    }

    pub mod shapes {
        use super::{grid_t_create, grid_t_unwrap, GridT};
        use crate::shapes;
        use crate::types::RealScalar;

        unsafe fn shapes_regular_sphere<T: RealScalar>(level: u32) -> *mut GridT {
            let wrapper = grid_t_create();
            let inner = unsafe { grid_t_unwrap(wrapper) }.unwrap();
            *inner = Box::new(shapes::regular_sphere::<T>(level));
            wrapper
        }
        #[no_mangle]
        pub unsafe extern "C" fn shapes_regular_sphere_f64(level: u32) -> *mut GridT {
            shapes_regular_sphere::<f64>(level)
        }
        #[no_mangle]
        pub unsafe extern "C" fn shapes_regular_sphere_f32(level: u32) -> *mut GridT {
            shapes_regular_sphere::<f32>(level)
        }
    }
}

pub mod entity {
    use super::EntityT;
    use super::{
        geometry_t_create, geometry_t_unwrap, topology_t_create, topology_t_unwrap, GeometryT,
        TopologyT,
    };
    use crate::{
        grid::serial::SingleElementGridEntity,
        traits::Entity,
        types::{Ownership, RealScalar},
    };
    use c_api_tools::{concretise_types, DType, DTypeIdentifier};
    use ndelement::{ciarlet::CiarletElement, types::ReferenceCellType};

    #[concretise_types(
        gen_type(name = "dtype", replace_with = ["f32", "f64"]),
        field(arg = 0, name = "wrap", wrapper = "EntityT", replace_with = ["SingleElementGridEntity<{{dtype}}, CiarletElement<{{dtype}}>>"]),
    )]
    pub fn entity_local_index<E: Entity>(entity: &E) -> usize {
        entity.local_index()
    }

    #[concretise_types(
        gen_type(name = "dtype", replace_with = ["f32", "f64"]),
        field(arg = 0, name = "wrap", wrapper = "EntityT", replace_with = ["SingleElementGridEntity<{{dtype}}, CiarletElement<{{dtype}}>>"]),
    )]
    pub fn entity_global_index<E: Entity>(entity: &E) -> usize {
        entity.global_index()
    }

    #[concretise_types(
        gen_type(name = "dtype", replace_with = ["f32", "f64"]),
        field(arg = 0, name = "wrap", wrapper = "EntityT", replace_with = ["SingleElementGridEntity<{{dtype}}, CiarletElement<{{dtype}}>>"]),
    )]
    pub fn entity_entity_type<E: Entity<EntityDescriptor = ReferenceCellType>>(
        entity: &E,
    ) -> ReferenceCellType {
        entity.entity_type()
    }

    #[concretise_types(
        gen_type(name = "dtype", replace_with = ["f32", "f64"]),
        field(arg = 0, name = "wrap", wrapper = "EntityT", replace_with = ["SingleElementGridEntity<{{dtype}}, CiarletElement<{{dtype}}>>"]),
    )]
    pub fn entity_has_id<E: Entity>(entity: &E) -> bool {
        entity.id().is_some()
    }

    #[concretise_types(
        gen_type(name = "dtype", replace_with = ["f32", "f64"]),
        field(arg = 0, name = "wrap", wrapper = "EntityT", replace_with = ["SingleElementGridEntity<{{dtype}}, CiarletElement<{{dtype}}>>"]),
    )]
    pub fn entity_id<E: Entity>(entity: &E) -> usize {
        entity.id().unwrap()
    }

    #[concretise_types(
        gen_type(name = "dtype", replace_with = ["f32", "f64"]),
        field(arg = 0, name = "wrap", wrapper = "EntityT", replace_with = ["SingleElementGridEntity<{{dtype}}, CiarletElement<{{dtype}}>>"]),
    )]
    pub fn entity_is_owned<E: Entity>(entity: &E) -> bool {
        entity.ownership() == Ownership::Owned
    }

    #[concretise_types(
        gen_type(name = "dtype", replace_with = ["f32", "f64"]),
        field(arg = 0, name = "wrap", wrapper = "EntityT", replace_with = ["SingleElementGridEntity<{{dtype}}, CiarletElement<{{dtype}}>>"]),
    )]
    pub fn entity_ownership_process<E: Entity>(entity: &E) -> usize {
        if let Ownership::Ghost(process, _index) = entity.ownership() {
            process
        } else {
            panic!("Cannot get owner of non-ghost entity");
        }
    }

    #[concretise_types(
        gen_type(name = "dtype", replace_with = ["f32", "f64"]),
        field(arg = 0, name = "wrap", wrapper = "EntityT", replace_with = ["SingleElementGridEntity<{{dtype}}, CiarletElement<{{dtype}}>>"]),
    )]
    pub fn entity_ownership_index<E: Entity>(entity: &E) -> usize {
        if let Ownership::Ghost(_process, index) = entity.ownership() {
            index
        } else {
            panic!("Cannot get owner of non-ghost entity");
        }
    }

    #[concretise_types(
        gen_type(name = "dtype", replace_with = ["f32", "f64"]),
        field(arg = 0, name = "wrap", wrapper = "EntityT", replace_with = ["SingleElementGridEntity<{{dtype}}, CiarletElement<{{dtype}}>>"]),
    )]
    pub fn entity_topology<E: Entity>(entity: &'static E) -> *mut TopologyT {
        let wrapper = topology_t_create();
        let inner = unsafe { topology_t_unwrap(wrapper) }.unwrap();
        *inner = Box::new(entity.topology());
        wrapper
    }

    #[concretise_types(
        gen_type(name = "dtype", replace_with = ["f32", "f64"]),
        field(arg = 0, name = "wrap", wrapper = "EntityT", replace_with = ["SingleElementGridEntity<{{dtype}}, CiarletElement<{{dtype}}>>"]),
    )]
    pub fn entity_geometry<E: Entity>(entity: &'static E) -> *mut GeometryT {
        let wrapper = geometry_t_create();
        let inner = unsafe { geometry_t_unwrap(wrapper) }.unwrap();
        *inner = Box::new(entity.geometry());
        wrapper
    }

    #[concretise_types(
        gen_type(name = "dtype", replace_with = ["f32", "f64"]),
        field(arg = 0, name = "wrap", wrapper = "EntityT", replace_with = ["SingleElementGridEntity<{{dtype}}, CiarletElement<{{dtype}}>>"]),
    )]
    pub fn entity_dtype<T: RealScalar + DTypeIdentifier, E: Entity<T = T>>(_entity: &E) -> DType {
        <T as DTypeIdentifier>::dtype()
    }
}

pub mod topology {
    use super::TopologyT;
    use crate::{topology::serial::SingleTypeEntityTopology, traits::Topology};
    use c_api_tools::concretise_types;

    #[concretise_types(
        gen_type(name = "", replace_with = [""]),
        field(arg = 0, name = "wrap", wrapper = "TopologyT", replace_with = ["SingleTypeEntityTopology"]),
    )]
    pub fn topology_sub_entity<T: Topology>(topology: &T, dim: usize, index: usize) -> usize {
        topology.sub_entity(dim, index)
    }

    #[concretise_types(
        gen_type(name = "", replace_with = [""]),
        field(arg = 0, name = "wrap", wrapper = "TopologyT", replace_with = ["SingleTypeEntityTopology"]),
    )]
    pub fn topology_sub_entities_size<T: Topology>(topology: &T, dim: usize) -> usize {
        topology.sub_entity_iter(dim).map(|_| 1).sum()
    }

    #[concretise_types(
        gen_type(name = "", replace_with = [""]),
        field(arg = 0, name = "wrap", wrapper = "TopologyT", replace_with = ["SingleTypeEntityTopology"]),
    )]
    pub fn topology_sub_entities<T: Topology>(topology: &T, dim: usize, entities: *mut usize) {
        for (i, e) in topology.sub_entity_iter(dim).enumerate() {
            unsafe {
                *entities.add(i) = e;
            }
        }
    }

    #[concretise_types(
        gen_type(name = "", replace_with = [""]),
        field(arg = 0, name = "wrap", wrapper = "TopologyT", replace_with = ["SingleTypeEntityTopology"]),
    )]
    pub fn topology_connected_entities_size<T: Topology>(topology: &T, dim: usize) -> usize {
        topology.connected_entity_iter(dim).map(|_| 1).sum()
    }

    #[concretise_types(
        gen_type(name = "", replace_with = [""]),
        field(arg = 0, name = "wrap", wrapper = "TopologyT", replace_with = ["SingleTypeEntityTopology"]),
    )]
    pub fn topology_connected_entities<T: Topology>(
        topology: &T,
        dim: usize,
        entities: *mut usize,
    ) {
        for (i, e) in topology.connected_entity_iter(dim).enumerate() {
            unsafe {
                *entities.add(i) = e;
            }
        }
    }
}

pub mod geometry {
    use super::GeometryT;
    use crate::{
        geometry::SingleElementEntityGeometry,
        traits::{Geometry, Point},
        types::RealScalar,
    };
    use c_api_tools::{concretise_types, DType, DTypeIdentifier};
    use ndelement::ciarlet::CiarletElement;
    use std::ffi::c_void;
    use std::slice::from_raw_parts_mut;

    #[concretise_types(
        gen_type(name = "dtype", replace_with = ["f32", "f64"]),
        field(arg = 0, name = "wrap", wrapper = "GeometryT", replace_with = ["SingleElementEntityGeometry<{{dtype}}, CiarletElement<{{dtype}}>>"]),
    )]
    pub fn geometry_points<G: Geometry>(geometry: &G, points: *mut c_void) {
        let points = points as *mut G::T;
        for (i, pt) in geometry.points().enumerate() {
            let gdim = pt.dim();
            pt.coords(unsafe { from_raw_parts_mut(points.add(gdim * i), gdim) });
        }
    }

    #[concretise_types(
        gen_type(name = "dtype", replace_with = ["f32", "f64"]),
        field(arg = 0, name = "wrap", wrapper = "GeometryT", replace_with = ["SingleElementEntityGeometry<{{dtype}}, CiarletElement<{{dtype}}>>"]),
    )]
    pub fn geometry_point_count<G: Geometry>(geometry: &G) -> usize {
        geometry.point_count()
    }

    #[concretise_types(
        gen_type(name = "dtype", replace_with = ["f32", "f64"]),
        field(arg = 0, name = "wrap", wrapper = "GeometryT", replace_with = ["SingleElementEntityGeometry<{{dtype}}, CiarletElement<{{dtype}}>>"]),
    )]
    pub fn geometry_degree<G: Geometry>(geometry: &G) -> usize {
        geometry.degree()
    }

    #[concretise_types(
        gen_type(name = "dtype", replace_with = ["f32", "f64"]),
        field(arg = 0, name = "wrap", wrapper = "GeometryT", replace_with = ["SingleElementEntityGeometry<{{dtype}}, CiarletElement<{{dtype}}>>"]),
    )]
    pub fn geometry_dtype<T: RealScalar + DTypeIdentifier, G: Geometry<T = T>>(
        _geometry: &G,
    ) -> DType {
        <T as DTypeIdentifier>::dtype()
    }
}

pub mod geometry_map {
    use super::GeometryMapT;
    use crate::{
        geometry::GeometryMap, traits::GeometryMap as GeometryMapTrait, types::RealScalar,
    };
    use c_api_tools::{concretise_types, DType, DTypeIdentifier};
    use std::ffi::c_void;
    use std::slice::from_raw_parts_mut;

    #[concretise_types(
        gen_type(name = "dtype", replace_with = ["f32", "f64"]),
        field(arg = 0, name = "wrap", wrapper = "GeometryMapT", replace_with = ["GeometryMap<{{dtype}}>"]),
    )]
    pub fn geometry_map_entity_topology_dimension<GM: GeometryMapTrait>(gmap: &GM) -> usize {
        gmap.entity_topology_dimension()
    }

    #[concretise_types(
        gen_type(name = "dtype", replace_with = ["f32", "f64"]),
        field(arg = 0, name = "wrap", wrapper = "GeometryMapT", replace_with = ["GeometryMap<{{dtype}}>"]),
    )]
    pub fn geometry_map_geometry_dimension<GM: GeometryMapTrait>(gmap: &GM) -> usize {
        gmap.geometry_dimension()
    }

    #[concretise_types(
        gen_type(name = "dtype", replace_with = ["f32", "f64"]),
        field(arg = 0, name = "wrap", wrapper = "GeometryMapT", replace_with = ["GeometryMap<{{dtype}}>"]),
    )]
    pub fn geometry_map_point_count<GM: GeometryMapTrait>(gmap: &GM) -> usize {
        gmap.point_count()
    }

    #[concretise_types(
        gen_type(name = "dtype", replace_with = ["f32", "f64"]),
        field(arg = 0, name = "wrap", wrapper = "GeometryMapT", replace_with = ["GeometryMap<{{dtype}}>"]),
    )]
    pub fn geometry_map_points<GM: GeometryMapTrait>(
        gmap: &GM,
        entity_index: usize,
        points: *mut c_void,
    ) {
        let points = points as *mut GM::T;
        gmap.points(entity_index, unsafe {
            from_raw_parts_mut(points, gmap.geometry_dimension() * gmap.point_count())
        });
    }

    #[concretise_types(
        gen_type(name = "dtype", replace_with = ["f32", "f64"]),
        field(arg = 0, name = "wrap", wrapper = "GeometryMapT", replace_with = ["GeometryMap<{{dtype}}>"]),
    )]
    pub fn geometry_map_jacobians<GM: GeometryMapTrait>(
        gmap: &GM,
        entity_index: usize,
        jacobians: *mut c_void,
    ) {
        let jacobians = jacobians as *mut GM::T;
        gmap.jacobians(entity_index, unsafe {
            from_raw_parts_mut(
                jacobians,
                gmap.geometry_dimension() * gmap.entity_topology_dimension() * gmap.point_count(),
            )
        });
    }
    #[concretise_types(
        gen_type(name = "dtype", replace_with = ["f32", "f64"]),
        field(arg = 0, name = "wrap", wrapper = "GeometryMapT", replace_with = ["GeometryMap<{{dtype}}>"]),
    )]
    pub fn geometry_map_jacobians_dets_normals<GM: GeometryMapTrait>(
        gmap: &GM,
        entity_index: usize,
        jacobians: *mut c_void,
        jdets: *mut c_void,
        normals: *mut c_void,
    ) {
        let jacobians = jacobians as *mut GM::T;
        let jdets = jdets as *mut GM::T;
        let normals = normals as *mut GM::T;
        gmap.jacobians_dets_normals(
            entity_index,
            unsafe {
                from_raw_parts_mut(
                    jacobians,
                    gmap.geometry_dimension()
                        * gmap.entity_topology_dimension()
                        * gmap.point_count(),
                )
            },
            unsafe { from_raw_parts_mut(jdets, gmap.point_count()) },
            unsafe { from_raw_parts_mut(normals, gmap.geometry_dimension() * gmap.point_count()) },
        );
    }

    #[concretise_types(
        gen_type(name = "dtype", replace_with = ["f32", "f64"]),
        field(arg = 0, name = "wrap", wrapper = "GeometryMapT", replace_with = ["GeometryMap<{{dtype}}>"]),
    )]
    pub fn geometry_map_dtype<T: RealScalar + DTypeIdentifier, GM: GeometryMapTrait<T = T>>(
        _gmap: &GM,
    ) -> DType {
        <T as DTypeIdentifier>::dtype()
    }
}
