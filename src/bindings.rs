//! Binding for C
#![allow(missing_docs)]
#![allow(clippy::missing_safety_doc)]

pub mod grid {
    use crate::{grid::SingleElementGrid, traits::Grid, types::RealScalar};
    use c_api_tools::DType;
    use c_api_tools::{cfuncs, concretise_types};
    use ndelement::{
        ciarlet::CiarletElement,
        ciarlet::LagrangeElementFamily,
        reference_cell,
        traits::{ElementFamily, FiniteElement},
        types::{Continuity, ReferenceCellType},
    };
    use std::ffi::c_void;

    /// Wrapper for grid types.
    #[cfuncs(name = "grid_t", create, free, unwrap)]
    pub struct GridT;

    /// Wrapper for entities.
    #[cfuncs(name = "entity_t", create, free, unwrap)]
    pub struct EntityT;

    #[allow(clippy::too_many_arguments)]
    pub unsafe fn single_element_grid_new(
        coordinates: *const c_void,
        npts: usize,
        gdim: usize,
        cells: *const usize,
        ncells: usize,
        cell_type: u8,
        geometry_degree: usize,
        dtype: DType,
    ) -> *mut GridT {
        pub fn single_elemeng_grid_impl<T: RealScalar>(
            coordinates: *const T,
            npts: usize,
            gdim: usize,
            cells: *const usize,
            ncells: usize,
            cell_type: u8,
            geometry_degree: usize,
        ) -> *mut GridT {
            let wrapper = grid_t_create();
            let inner = unsafe { grid_t_unwrap(wrapper) }.unwrap();
            let coordinates = unsafe { std::slice::from_raw_parts(coordinates, npts * gdim) };
            let points_per_cell =
                LagrangeElementFamily::<T>::new(geometry_degree, Continuity::Standard)
                    .element(ReferenceCellType::from(cell_type).unwrap())
                    .dim();
            let cells = unsafe { std::slice::from_raw_parts(cells, ncells * points_per_cell) };

            *inner = Box::new(
                SingleElementGrid::<T, CiarletElement<T>>::new_from_raw_data(
                    coordinates,
                    gdim,
                    cells,
                    ReferenceCellType::from(cell_type).unwrap(),
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
    pub fn grid_tdim<T: RealScalar>(grid: &SingleElementGrid<T, CiarletElement<T>>) -> usize {
        grid.topology_dim()
    }

    #[concretise_types(
        gen_type(name = "dtype", replace_with = ["f32", "f64"]),
        field(arg = 0, name = "wrap", wrapper = "GridT", replace_with = ["SingleElementGrid<{{dtype}}, CiarletElement<{{dtype}}>>"]),
    )]
    pub fn grid_gdim<T: RealScalar>(grid: &SingleElementGrid<T, CiarletElement<T>>) -> usize {
        grid.geometry_dim()
    }

    #[concretise_types(
        gen_type(name = "dtype", replace_with = ["f32", "f64"]),
        field(arg = 0, name = "wrap", wrapper = "GridT", replace_with = ["SingleElementGrid<{{dtype}}, CiarletElement<{{dtype}}>>"]),
    )]
    pub fn grid_entity_count<T: RealScalar>(
        grid: &SingleElementGrid<T, CiarletElement<T>>,
        entity_type: ReferenceCellType,
    ) -> usize {
        grid.entity_count(entity_type)
    }

    #[concretise_types(
        gen_type(name = "dtype", replace_with = ["f32", "f64"]),
        field(arg = 0, name = "wrap", wrapper = "GridT", replace_with = ["SingleElementGrid<{{dtype}}, CiarletElement<{{dtype}}>>"]),
    )]
    pub fn grid_entity<T: RealScalar>(
        grid: &'static SingleElementGrid<T, CiarletElement<T>>,
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

    pub fn grid_entity_from_id<T: RealScalar>(
        grid: &'static SingleElementGrid<T, CiarletElement<T>>,
        dim: usize,
        id: usize,
    ) -> *const EntityT {
        let wrapper = entity_t_create();
        let inner = unsafe { entity_t_unwrap(wrapper) }.unwrap();
        let entity = grid.entity_from_id(dim, id);

        if let Some(entity) = entity {
            *inner = Box::new(entity);
        }

        wrapper
    }
}

// pub mod grid {
//     use super::entity::{EntityType, EntityWrapper};
//     use super::geometry_map::GeometryMapWrapper;
//     use super::DType;
//     use crate::{grid::SingleElementGrid, traits::Grid, types::RealScalar};
//     use ndelement::{
//         ciarlet::CiarletElement,
//         ciarlet::LagrangeElementFamily,
//         reference_cell,
//         traits::{ElementFamily, FiniteElement},
//         types::{Continuity, ReferenceCellType},
//     };
//     use std::ffi::c_void;
//     use std::slice::from_raw_parts;

//     unsafe fn grid_entity_from_id_internal<T: Grid>(
//         grid: *const GridWrapper,
//         dim: usize,
//         id: usize,
//         etype: EntityType,
//     ) -> *const EntityWrapper {
//         let entity = EntityWrapper {
//             entity: Box::into_raw(Box::new(
//                 (*extract_grid::<T>(grid)).entity_from_id(dim, id).unwrap(),
//             )) as *const c_void,
//             dtype: (*grid).dtype,
//             etype,
//         };
//         Box::into_raw(Box::new(entity))
//     }

//     unsafe fn grid_geometry_map_internal<
//         T: RealScalar,
//         G: Grid<EntityDescriptor = ReferenceCellType, T = T>,
//     >(
//         grid: *const GridWrapper,
//         entity_type: u8,
//         points: *const T,
//         npoints: usize,
//     ) -> *const GeometryMapWrapper {
//         let etype = ReferenceCellType::from(entity_type).unwrap();
//         let gmap = GeometryMapWrapper {
//             geometry_map: Box::into_raw(Box::new((*extract_grid::<G>(grid)).geometry_map(
//                 etype,
//                 from_raw_parts(points, npoints * reference_cell::dim(etype)),
//             ))) as *const c_void,
//             dtype: (*grid).dtype,
//         };
//         Box::into_raw(Box::new(gmap))
//     }

//     #[no_mangle]
//     pub unsafe extern "C" fn grid_geometry_map(
//         grid: *const GridWrapper,
//         entity_type: u8,
//         points: *const c_void,
//         npoints: usize,
//     ) -> *const GeometryMapWrapper {
//         match (*grid).gtype {
//             GridType::SerialSingleElementGrid => match (*grid).dtype {
//                 DType::F32 => grid_geometry_map_internal::<
//                     f32,
//                     SingleElementGrid<f32, CiarletElement<f32>>,
//                 >(grid, entity_type, points as *const f32, npoints),
//                 DType::F64 => grid_geometry_map_internal::<
//                     f64,
//                     SingleElementGrid<f64, CiarletElement<f64>>,
//                 >(grid, entity_type, points as *const f64, npoints),
//             },
//         }
//     }

//     #[no_mangle]
//     pub unsafe extern "C" fn grid_entity_types_size(grid: *const GridWrapper, dim: usize) -> usize {
//         match (*grid).gtype {
//             GridType::SerialSingleElementGrid => match (*grid).dtype {
//                 DType::F32 => (*extract_grid::<SingleElementGrid<f32, CiarletElement<f32>>>(grid))
//                     .entity_types(dim),
//                 DType::F64 => (*extract_grid::<SingleElementGrid<f64, CiarletElement<f64>>>(grid))
//                     .entity_types(dim),
//             },
//         }
//         .len()
//     }

//     #[no_mangle]
//     pub unsafe extern "C" fn grid_entity_types(
//         grid: *const GridWrapper,
//         dim: usize,
//         types: *mut u8,
//     ) {
//         for (i, t) in match (*grid).gtype {
//             GridType::SerialSingleElementGrid => match (*grid).dtype {
//                 DType::F32 => (*extract_grid::<SingleElementGrid<f32, CiarletElement<f32>>>(grid))
//                     .entity_types(dim),
//                 DType::F64 => (*extract_grid::<SingleElementGrid<f64, CiarletElement<f64>>>(grid))
//                     .entity_types(dim),
//             },
//         }
//         .iter()
//         .enumerate()
//         {
//             *types.add(i) = *t as u8;
//         }
//     }

//     #[no_mangle]
//     pub unsafe extern "C" fn grid_dtype(grid: *const GridWrapper) -> u8 {
//         (*grid).dtype as u8
//     }
// }

// pub mod entity {
//     use super::geometry::{GeometryType, GeometryWrapper};
//     use super::topology::{TopologyType, TopologyWrapper};
//     use super::DType;
//     use crate::{grid::serial::SingleElementGridEntity, traits::Entity, types::Ownership};
//     use ndelement::ciarlet::CiarletElement;
//     use std::ffi::c_void;

//     #[derive(Debug, PartialEq, Clone, Copy)]
//     #[repr(u8)]
//     pub enum EntityType {
//         SingleElementGridEntity = 0,
//     }

//     #[repr(C)]
//     pub struct EntityWrapper {
//         pub entity: *const c_void,
//         pub dtype: DType,
//         pub etype: EntityType,
//     }

//     impl Drop for EntityWrapper {
//         fn drop(&mut self) {
//             let Self {
//                 entity,
//                 dtype,
//                 etype,
//             } = self;
//             match etype {
//                 EntityType::SingleElementGridEntity => match dtype {
//                     DType::F32 => drop(unsafe {
//                         Box::from_raw(
//                             *entity as *mut SingleElementGridEntity<f32, CiarletElement<f32>>,
//                         )
//                     }),
//                     DType::F64 => drop(unsafe {
//                         Box::from_raw(
//                             *entity as *mut SingleElementGridEntity<f64, CiarletElement<f64>>,
//                         )
//                     }),
//                 },
//             }
//         }
//     }

//     #[no_mangle]
//     pub unsafe extern "C" fn free_entity(e: *mut EntityWrapper) {
//         assert!(!e.is_null());
//         unsafe { drop(Box::from_raw(e)) }
//     }

//     unsafe fn extract_entity<E: Entity>(entity: *const EntityWrapper) -> *const E {
//         (*entity).entity as *const E
//     }

//     #[no_mangle]
//     pub unsafe extern "C" fn entity_local_index(entity: *const EntityWrapper) -> usize {
//         match (*entity).etype {
//             EntityType::SingleElementGridEntity => match (*entity).dtype {
//                 DType::F32 => {
//                     (*extract_entity::<SingleElementGridEntity<f32, CiarletElement<f32>>>(entity))
//                         .local_index()
//                 }
//                 DType::F64 => {
//                     (*extract_entity::<SingleElementGridEntity<f64, CiarletElement<f64>>>(entity))
//                         .local_index()
//                 }
//             },
//         }
//     }

//     #[no_mangle]
//     pub unsafe extern "C" fn entity_global_index(entity: *const EntityWrapper) -> usize {
//         match (*entity).etype {
//             EntityType::SingleElementGridEntity => match (*entity).dtype {
//                 DType::F32 => {
//                     (*extract_entity::<SingleElementGridEntity<f32, CiarletElement<f32>>>(entity))
//                         .global_index()
//                 }
//                 DType::F64 => {
//                     (*extract_entity::<SingleElementGridEntity<f64, CiarletElement<f64>>>(entity))
//                         .global_index()
//                 }
//             },
//         }
//     }

//     #[no_mangle]
//     pub unsafe extern "C" fn entity_entity_type(entity: *const EntityWrapper) -> u8 {
//         match (*entity).etype {
//             EntityType::SingleElementGridEntity => match (*entity).dtype {
//                 DType::F32 => {
//                     (*extract_entity::<SingleElementGridEntity<f32, CiarletElement<f32>>>(entity))
//                         .entity_type() as u8
//                 }

//                 DType::F64 => {
//                     (*extract_entity::<SingleElementGridEntity<f64, CiarletElement<f64>>>(entity))
//                         .entity_type() as u8
//                 }
//             },
//         }
//     }
//     #[no_mangle]
//     pub unsafe extern "C" fn entity_has_id(entity: *const EntityWrapper) -> bool {
//         match (*entity).etype {
//             EntityType::SingleElementGridEntity => match (*entity).dtype {
//                 DType::F32 => {
//                     (*extract_entity::<SingleElementGridEntity<f32, CiarletElement<f32>>>(entity))
//                         .id()
//                 }
//                 DType::F64 => {
//                     (*extract_entity::<SingleElementGridEntity<f64, CiarletElement<f64>>>(entity))
//                         .id()
//                 }
//             },
//         }
//         .is_some()
//     }

//     #[no_mangle]
//     pub unsafe extern "C" fn entity_id(entity: *const EntityWrapper) -> usize {
//         match (*entity).etype {
//             EntityType::SingleElementGridEntity => match (*entity).dtype {
//                 DType::F32 => {
//                     (*extract_entity::<SingleElementGridEntity<f32, CiarletElement<f32>>>(entity))
//                         .id()
//                 }
//                 DType::F64 => {
//                     (*extract_entity::<SingleElementGridEntity<f64, CiarletElement<f64>>>(entity))
//                         .id()
//                 }
//             },
//         }
//         .unwrap()
//     }

//     #[no_mangle]
//     pub unsafe extern "C" fn entity_is_owned(entity: *const EntityWrapper) -> bool {
//         let ownership = match (*entity).etype {
//             EntityType::SingleElementGridEntity => match (*entity).dtype {
//                 DType::F32 => {
//                     (*extract_entity::<SingleElementGridEntity<f32, CiarletElement<f32>>>(entity))
//                         .ownership()
//                 }
//                 DType::F64 => {
//                     (*extract_entity::<SingleElementGridEntity<f64, CiarletElement<f64>>>(entity))
//                         .ownership()
//                 }
//             },
//         };

//         ownership == Ownership::Owned
//     }

//     #[no_mangle]
//     pub unsafe extern "C" fn entity_ownership_process(entity: *const EntityWrapper) -> usize {
//         let ownership = match (*entity).etype {
//             EntityType::SingleElementGridEntity => match (*entity).dtype {
//                 DType::F32 => {
//                     (*extract_entity::<SingleElementGridEntity<f32, CiarletElement<f32>>>(entity))
//                         .ownership()
//                 }
//                 DType::F64 => {
//                     (*extract_entity::<SingleElementGridEntity<f64, CiarletElement<f64>>>(entity))
//                         .ownership()
//                 }
//             },
//         };

//         if let Ownership::Ghost(process, _index) = ownership {
//             process
//         } else {
//             panic!();
//         }
//     }

//     #[no_mangle]
//     pub unsafe extern "C" fn entity_ownership_index(entity: *const EntityWrapper) -> usize {
//         let ownership = match (*entity).etype {
//             EntityType::SingleElementGridEntity => match (*entity).dtype {
//                 DType::F32 => {
//                     (*extract_entity::<SingleElementGridEntity<f32, CiarletElement<f32>>>(entity))
//                         .ownership()
//                 }
//                 DType::F64 => {
//                     (*extract_entity::<SingleElementGridEntity<f64, CiarletElement<f64>>>(entity))
//                         .ownership()
//                 }
//             },
//         };

//         if let Ownership::Ghost(_process, index) = ownership {
//             index
//         } else {
//             panic!();
//         }
//     }

//     unsafe fn entity_topology_internal<E: Entity>(
//         entity: *const EntityWrapper,
//         ttype: TopologyType,
//     ) -> *const TopologyWrapper {
//         let topology = TopologyWrapper {
//             topology: Box::into_raw(Box::new((*extract_entity::<E>(entity)).topology()))
//                 as *const c_void,
//             ttype,
//         };
//         Box::into_raw(Box::new(topology))
//     }

//     #[no_mangle]
//     pub unsafe extern "C" fn entity_topology(
//         entity: *const EntityWrapper,
//     ) -> *const TopologyWrapper {
//         match (*entity).etype {
//             EntityType::SingleElementGridEntity => match (*entity).dtype {
//                 DType::F32 => entity_topology_internal::<
//                     SingleElementGridEntity<f32, CiarletElement<f32>>,
//                 >(entity, TopologyType::SingleTypeEntityTopology),
//                 DType::F64 => entity_topology_internal::<
//                     SingleElementGridEntity<f64, CiarletElement<f64>>,
//                 >(entity, TopologyType::SingleTypeEntityTopology),
//             },
//         }
//     }

//     unsafe fn entity_geometry_internal<E: Entity>(
//         entity: *const EntityWrapper,
//         gtype: GeometryType,
//         dtype: DType,
//     ) -> *const GeometryWrapper {
//         let geometry = GeometryWrapper {
//             geometry: Box::into_raw(Box::new((*extract_entity::<E>(entity)).geometry()))
//                 as *const c_void,
//             gtype,
//             dtype,
//         };
//         Box::into_raw(Box::new(geometry))
//     }

//     #[no_mangle]
//     pub unsafe extern "C" fn entity_geometry(
//         entity: *const EntityWrapper,
//     ) -> *const GeometryWrapper {
//         match (*entity).etype {
//             EntityType::SingleElementGridEntity => match (*entity).dtype {
//                 DType::F32 => {
//                     entity_geometry_internal::<SingleElementGridEntity<f32, CiarletElement<f32>>>(
//                         entity,
//                         GeometryType::SingleElementEntityGeometry,
//                         (*entity).dtype,
//                     )
//                 }
//                 DType::F64 => {
//                     entity_geometry_internal::<SingleElementGridEntity<f64, CiarletElement<f64>>>(
//                         entity,
//                         GeometryType::SingleElementEntityGeometry,
//                         (*entity).dtype,
//                     )
//                 }
//             },
//         }
//     }
//     #[no_mangle]
//     pub unsafe extern "C" fn entity_dtype(entity: *const EntityWrapper) -> u8 {
//         (*entity).dtype as u8
//     }
// }

// pub mod topology {
//     use crate::{topology::serial::SingleTypeEntityTopology, traits::Topology};
//     use std::ffi::c_void;

//     #[derive(Debug, PartialEq, Clone, Copy)]
//     #[repr(u8)]
//     pub enum TopologyType {
//         SingleTypeEntityTopology = 0,
//     }

//     #[repr(C)]
//     pub struct TopologyWrapper {
//         pub topology: *const c_void,
//         pub ttype: TopologyType,
//     }

//     impl Drop for TopologyWrapper {
//         fn drop(&mut self) {
//             let Self { topology, ttype } = self;
//             match ttype {
//                 TopologyType::SingleTypeEntityTopology => {
//                     drop(unsafe { Box::from_raw(*topology as *mut SingleTypeEntityTopology) })
//                 }
//             }
//         }
//     }

//     #[no_mangle]
//     pub unsafe extern "C" fn free_topology(t: *mut TopologyWrapper) {
//         assert!(!t.is_null());
//         unsafe { drop(Box::from_raw(t)) }
//     }

//     unsafe fn extract_topology<T: Topology>(topology: *const TopologyWrapper) -> *const T {
//         (*topology).topology as *const T
//     }

//     #[no_mangle]
//     pub unsafe extern "C" fn topology_sub_entity(
//         topology: *const TopologyWrapper,
//         dim: usize,
//         index: usize,
//     ) -> usize {
//         match (*topology).ttype {
//             TopologyType::SingleTypeEntityTopology => {
//                 (*extract_topology::<SingleTypeEntityTopology>(topology)).sub_entity(dim, index)
//             }
//         }
//     }

//     #[no_mangle]
//     pub unsafe extern "C" fn topology_sub_entities_size(
//         topology: *const TopologyWrapper,
//         dim: usize,
//     ) -> usize {
//         match (*topology).ttype {
//             TopologyType::SingleTypeEntityTopology => {
//                 (*extract_topology::<SingleTypeEntityTopology>(topology)).sub_entity_iter(dim)
//             }
//         }
//         .len()
//     }

//     #[no_mangle]
//     pub unsafe extern "C" fn topology_sub_entities(
//         topology: *const TopologyWrapper,
//         dim: usize,
//         entities: *mut usize,
//     ) {
//         for (i, e) in match (*topology).ttype {
//             TopologyType::SingleTypeEntityTopology => {
//                 (*extract_topology::<SingleTypeEntityTopology>(topology)).sub_entity_iter(dim)
//             }
//         }
//         .enumerate()
//         {
//             *entities.add(i) = e;
//         }
//     }

//     #[no_mangle]
//     pub unsafe extern "C" fn topology_connected_entities_size(
//         topology: *const TopologyWrapper,
//         dim: usize,
//     ) -> usize {
//         match (*topology).ttype {
//             TopologyType::SingleTypeEntityTopology => {
//                 (*extract_topology::<SingleTypeEntityTopology>(topology)).connected_entity_iter(dim)
//             }
//         }
//         .len()
//     }

//     #[no_mangle]
//     pub unsafe extern "C" fn topology_connected_entities(
//         topology: *const TopologyWrapper,
//         dim: usize,
//         entities: *mut usize,
//     ) {
//         for (i, e) in match (*topology).ttype {
//             TopologyType::SingleTypeEntityTopology => {
//                 (*extract_topology::<SingleTypeEntityTopology>(topology)).connected_entity_iter(dim)
//             }
//         }
//         .enumerate()
//         {
//             *entities.add(i) = e;
//         }
//     }
// }

// pub mod geometry {
//     use super::DType;
//     use crate::{
//         geometry::SingleElementEntityGeometry,
//         traits::{Geometry, Point},
//         types::RealScalar,
//     };
//     use ndelement::ciarlet::CiarletElement;
//     use std::ffi::c_void;
//     use std::slice::from_raw_parts_mut;

//     #[derive(Debug, PartialEq, Clone, Copy)]
//     #[repr(u8)]
//     pub enum GeometryType {
//         SingleElementEntityGeometry = 0,
//     }

//     #[repr(C)]
//     pub struct GeometryWrapper {
//         pub geometry: *const c_void,
//         pub dtype: DType,
//         pub gtype: GeometryType,
//     }

//     impl Drop for GeometryWrapper {
//         fn drop(&mut self) {
//             let Self {
//                 geometry,
//                 dtype,
//                 gtype,
//             } = self;
//             match gtype {
//                 GeometryType::SingleElementEntityGeometry => match dtype {
//                     DType::F32 => drop(unsafe {
//                         Box::from_raw(
//                             *geometry as *mut SingleElementEntityGeometry<f32, CiarletElement<f32>>,
//                         )
//                     }),
//                     DType::F64 => drop(unsafe {
//                         Box::from_raw(
//                             *geometry as *mut SingleElementEntityGeometry<f64, CiarletElement<f64>>,
//                         )
//                     }),
//                 },
//             }
//         }
//     }

//     #[no_mangle]
//     pub unsafe extern "C" fn free_geometry(g: *mut GeometryWrapper) {
//         assert!(!g.is_null());
//         unsafe { drop(Box::from_raw(g)) }
//     }

//     unsafe fn extract_geometry<G: Geometry>(geometry: *const GeometryWrapper) -> *const G {
//         (*geometry).geometry as *const G
//     }

//     unsafe fn geometry_points_internal<T: RealScalar, G: Geometry<T = T>>(
//         geometry: *mut GeometryWrapper,
//         points: *mut c_void,
//     ) {
//         let points = points as *mut T;
//         for (i, pt) in (*extract_geometry::<G>(geometry)).points().enumerate() {
//             let gdim = pt.dim();
//             pt.coords(from_raw_parts_mut(points.add(gdim * i), gdim));
//         }
//     }

//     #[no_mangle]
//     pub unsafe extern "C" fn geometry_points(geometry: *mut GeometryWrapper, points: *mut c_void) {
//         match (*geometry).gtype {
//             GeometryType::SingleElementEntityGeometry => match (*geometry).dtype {
//                 DType::F32 => geometry_points_internal::<
//                     f32,
//                     SingleElementEntityGeometry<f32, CiarletElement<f32>>,
//                 >(geometry, points),
//                 DType::F64 => geometry_points_internal::<
//                     f64,
//                     SingleElementEntityGeometry<f64, CiarletElement<f64>>,
//                 >(geometry, points),
//             },
//         }
//     }

//     #[no_mangle]
//     pub unsafe extern "C" fn geometry_point_count(geometry: *mut GeometryWrapper) -> usize {
//         match (*geometry).gtype {
//             GeometryType::SingleElementEntityGeometry => match (*geometry).dtype {
//                 DType::F32 => (*extract_geometry::<
//                     SingleElementEntityGeometry<f32, CiarletElement<f32>>,
//                 >(geometry))
//                 .point_count(),
//                 DType::F64 => (*extract_geometry::<
//                     SingleElementEntityGeometry<f64, CiarletElement<f64>>,
//                 >(geometry))
//                 .point_count(),
//             },
//         }
//     }

//     #[no_mangle]
//     pub unsafe extern "C" fn geometry_degree(geometry: *mut GeometryWrapper) -> usize {
//         match (*geometry).gtype {
//             GeometryType::SingleElementEntityGeometry => match (*geometry).dtype {
//                 DType::F32 => (*extract_geometry::<
//                     SingleElementEntityGeometry<f32, CiarletElement<f32>>,
//                 >(geometry))
//                 .degree(),
//                 DType::F64 => (*extract_geometry::<
//                     SingleElementEntityGeometry<f64, CiarletElement<f64>>,
//                 >(geometry))
//                 .degree(),
//             },
//         }
//     }

//     #[no_mangle]
//     pub unsafe extern "C" fn geometry_dtype(geometry: *const GeometryWrapper) -> u8 {
//         (*geometry).dtype as u8
//     }
// }

// pub mod geometry_map {
//     use super::DType;
//     use crate::{
//         geometry::GeometryMap, traits::GeometryMap as GeometryMapTrait, types::RealScalar,
//     };
//     use std::ffi::c_void;
//     use std::slice::from_raw_parts_mut;

//     #[repr(C)]
//     pub struct GeometryMapWrapper {
//         pub geometry_map: *const c_void,
//         pub dtype: DType,
//     }

//     impl Drop for GeometryMapWrapper {
//         fn drop(&mut self) {
//             let Self {
//                 geometry_map,
//                 dtype,
//             } = self;
//             match dtype {
//                 DType::F32 => {
//                     drop(unsafe { Box::from_raw(*geometry_map as *mut GeometryMap<f32>) })
//                 }
//                 DType::F64 => {
//                     drop(unsafe { Box::from_raw(*geometry_map as *mut GeometryMap<f64>) })
//                 }
//             }
//         }
//     }

//     #[no_mangle]
//     pub unsafe extern "C" fn free_geometry_map(g: *mut GeometryMapWrapper) {
//         assert!(!g.is_null());
//         unsafe { drop(Box::from_raw(g)) }
//     }

//     unsafe fn extract_geometry_map<'a, T: RealScalar>(
//         gmap: *const GeometryMapWrapper,
//     ) -> *const GeometryMap<'a, T> {
//         (*gmap).geometry_map as *const GeometryMap<T>
//     }

//     #[no_mangle]
//     pub unsafe extern "C" fn geometry_map_entity_topology_dimension(
//         gmap: *mut GeometryMapWrapper,
//     ) -> usize {
//         match (*gmap).dtype {
//             DType::F32 => (*extract_geometry_map::<f32>(gmap)).entity_topology_dimension(),
//             DType::F64 => (*extract_geometry_map::<f64>(gmap)).entity_topology_dimension(),
//         }
//     }

//     #[no_mangle]
//     pub unsafe extern "C" fn geometry_map_geometry_dimension(
//         gmap: *mut GeometryMapWrapper,
//     ) -> usize {
//         match (*gmap).dtype {
//             DType::F32 => (*extract_geometry_map::<f32>(gmap)).geometry_dimension(),
//             DType::F64 => (*extract_geometry_map::<f64>(gmap)).geometry_dimension(),
//         }
//     }

//     #[no_mangle]
//     pub unsafe extern "C" fn geometry_map_point_count(gmap: *mut GeometryMapWrapper) -> usize {
//         match (*gmap).dtype {
//             DType::F32 => (*extract_geometry_map::<f32>(gmap)).point_count(),
//             DType::F64 => (*extract_geometry_map::<f64>(gmap)).point_count(),
//         }
//     }

//     unsafe fn geometry_map_points_internal<T: RealScalar>(
//         gmap: *mut GeometryMapWrapper,
//         entity_index: usize,
//         points: *mut T,
//     ) {
//         let map = extract_geometry_map::<T>(gmap);
//         (*map).points(
//             entity_index,
//             from_raw_parts_mut(points, (*map).geometry_dimension() * (*map).point_count()),
//         );
//     }

//     #[no_mangle]
//     pub unsafe extern "C" fn geometry_map_points(
//         gmap: *mut GeometryMapWrapper,
//         entity_index: usize,
//         points: *mut c_void,
//     ) {
//         match (*gmap).dtype {
//             DType::F32 => {
//                 geometry_map_points_internal::<f32>(gmap, entity_index, points as *mut f32)
//             }
//             DType::F64 => {
//                 geometry_map_points_internal::<f64>(gmap, entity_index, points as *mut f64)
//             }
//         }
//     }

//     unsafe fn geometry_map_jacobians_internal<T: RealScalar>(
//         gmap: *mut GeometryMapWrapper,
//         entity_index: usize,
//         jacobians: *mut T,
//     ) {
//         let map = extract_geometry_map::<T>(gmap);
//         (*map).jacobians(
//             entity_index,
//             from_raw_parts_mut(
//                 jacobians,
//                 (*map).geometry_dimension()
//                     * (*map).entity_topology_dimension()
//                     * (*map).point_count(),
//             ),
//         );
//     }

//     #[no_mangle]
//     pub unsafe extern "C" fn geometry_map_jacobians(
//         gmap: *mut GeometryMapWrapper,
//         entity_index: usize,
//         jacobians: *mut c_void,
//     ) {
//         match (*gmap).dtype {
//             DType::F32 => {
//                 geometry_map_jacobians_internal::<f32>(gmap, entity_index, jacobians as *mut f32)
//             }
//             DType::F64 => {
//                 geometry_map_jacobians_internal::<f64>(gmap, entity_index, jacobians as *mut f64)
//             }
//         }
//     }

//     unsafe fn geometry_map_jacobians_dets_normals_internal<T: RealScalar>(
//         gmap: *mut GeometryMapWrapper,
//         entity_index: usize,
//         jacobians: *mut T,
//         jdets: *mut T,
//         normals: *mut T,
//     ) {
//         let map = extract_geometry_map::<T>(gmap);
//         (*map).jacobians_dets_normals(
//             entity_index,
//             from_raw_parts_mut(
//                 jacobians,
//                 (*map).geometry_dimension()
//                     * (*map).entity_topology_dimension()
//                     * (*map).point_count(),
//             ),
//             from_raw_parts_mut(jdets, (*map).point_count()),
//             from_raw_parts_mut(normals, (*map).geometry_dimension() * (*map).point_count()),
//         );
//     }

//     #[no_mangle]
//     pub unsafe extern "C" fn geometry_map_jacobians_dets_normals(
//         gmap: *mut GeometryMapWrapper,
//         entity_index: usize,
//         jacobians: *mut c_void,
//         jdets: *mut c_void,
//         normals: *mut c_void,
//     ) {
//         match (*gmap).dtype {
//             DType::F32 => geometry_map_jacobians_dets_normals_internal::<f32>(
//                 gmap,
//                 entity_index,
//                 jacobians as *mut f32,
//                 jdets as *mut f32,
//                 normals as *mut f32,
//             ),
//             DType::F64 => geometry_map_jacobians_dets_normals_internal::<f64>(
//                 gmap,
//                 entity_index,
//                 jacobians as *mut f64,
//                 jdets as *mut f64,
//                 normals as *mut f64,
//             ),
//         }
//     }

//     #[no_mangle]
//     pub unsafe extern "C" fn geometry_map_dtype(gmap: *const GeometryMapWrapper) -> u8 {
//         (*gmap).dtype as u8
//     }
// }

// pub mod shapes {
//     use super::grid::{GridType, GridWrapper};
//     use super::DType;
//     use crate::shapes::regular_sphere;
//     use crate::types::RealScalar;
//     use std::ffi::c_void;

//     unsafe fn shapes_regular_sphere<T: RealScalar>(level: u32, dtype: DType) -> *const GridWrapper {
//         let grid = Box::into_raw(Box::new(regular_sphere::<T>(level))) as *const c_void;
//         Box::into_raw(Box::new(GridWrapper {
//             grid,
//             dtype,
//             gtype: GridType::SerialSingleElementGrid,
//         }))
//     }
//     #[no_mangle]
//     pub unsafe extern "C" fn shapes_regular_sphere_f64(level: u32) -> *const GridWrapper {
//         shapes_regular_sphere::<f64>(level, DType::F64)
//     }
//     #[no_mangle]
//     pub unsafe extern "C" fn shapes_regular_sphere_f32(level: u32) -> *const GridWrapper {
//         shapes_regular_sphere::<f32>(level, DType::F32)
//     }
// }
