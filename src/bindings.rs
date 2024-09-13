//! Binding for C
#![allow(missing_docs)]
#![allow(clippy::missing_safety_doc)]

#[derive(Debug, PartialEq, Clone, Copy)]
#[repr(u8)]
pub enum DType {
    F32 = 0,
    F64 = 1,
}

mod grid {
    use super::entity::{EntityType, EntityWrapper};
    use super::DType;
    use crate::{grid::SingleElementGrid, traits::Grid, types::RealScalar};
    use ndelement::{
        ciarlet::CiarletElement,
        ciarlet::LagrangeElementFamily,
        traits::{ElementFamily, FiniteElement},
        types::{Continuity, ReferenceCellType},
    };
    use std::ffi::c_void;
    use std::slice::from_raw_parts;

    #[derive(Debug, PartialEq, Clone, Copy)]
    #[repr(u8)]
    pub enum GridType {
        SerialSingleElementGrid = 0,
    }

    #[repr(C)]
    pub struct GridWrapper {
        pub grid: *const c_void,
        pub dtype: DType,
        pub gtype: GridType,
    }

    impl Drop for GridWrapper {
        fn drop(&mut self) {
            let Self { grid, dtype, gtype } = self;
            match gtype {
                GridType::SerialSingleElementGrid => match dtype {
                    DType::F32 => drop(unsafe {
                        Box::from_raw(*grid as *mut SingleElementGrid<f32, CiarletElement<f32>>)
                    }),
                    DType::F64 => drop(unsafe {
                        Box::from_raw(*grid as *mut SingleElementGrid<f64, CiarletElement<f64>>)
                    }),
                },
            }
        }
    }

    #[no_mangle]
    pub unsafe extern "C" fn free_grid(g: *mut GridWrapper) {
        assert!(!g.is_null());
        unsafe { drop(Box::from_raw(g)) }
    }

    #[allow(clippy::too_many_arguments)]
    unsafe fn single_element_grid_new_from_raw_data<T: RealScalar>(
        coordinates: *const T,
        npts: usize,
        gdim: usize,
        cells: *const usize,
        ncells: usize,
        cell_type: u8,
        geometry_degree: usize,
        dtype: DType,
    ) -> *const GridWrapper {
        let coordinates = from_raw_parts(coordinates, npts * gdim);
        let points_per_cell =
            LagrangeElementFamily::<T>::new(geometry_degree, Continuity::Standard)
                .element(ReferenceCellType::from(cell_type).unwrap())
                .dim();
        let cells = from_raw_parts(cells, ncells * points_per_cell);

        let grid = Box::into_raw(Box::new(
            SingleElementGrid::<T, CiarletElement<T>>::new_from_raw_data(
                coordinates,
                gdim,
                cells,
                ReferenceCellType::from(cell_type).unwrap(),
                geometry_degree,
            ),
        )) as *const c_void;
        Box::into_raw(Box::new(GridWrapper {
            grid,
            dtype,
            gtype: GridType::SerialSingleElementGrid,
        }))
    }

    #[no_mangle]
    pub unsafe extern "C" fn single_element_grid_new_from_raw_data_f32(
        coordinates: *const f32,
        npts: usize,
        gdim: usize,
        cells: *const usize,
        ncells: usize,
        cell_type: u8,
        geometry_degree: usize,
    ) -> *const GridWrapper {
        single_element_grid_new_from_raw_data::<f32>(
            coordinates,
            npts,
            gdim,
            cells,
            ncells,
            cell_type,
            geometry_degree,
            DType::F32,
        )
    }

    #[no_mangle]
    pub unsafe extern "C" fn single_element_grid_new_from_raw_data_f64(
        coordinates: *const f64,
        npts: usize,
        gdim: usize,
        cells: *const usize,
        ncells: usize,
        cell_type: u8,
        geometry_degree: usize,
    ) -> *const GridWrapper {
        single_element_grid_new_from_raw_data::<f64>(
            coordinates,
            npts,
            gdim,
            cells,
            ncells,
            cell_type,
            geometry_degree,
            DType::F64,
        )
    }

    unsafe fn extract_grid<G: Grid>(grid: *const GridWrapper) -> *const G {
        (*grid).grid as *const G
    }

    #[no_mangle]
    pub unsafe extern "C" fn grid_tdim(grid: *const GridWrapper) -> usize {
        match (*grid).gtype {
            GridType::SerialSingleElementGrid => match (*grid).dtype {
                DType::F32 => (*extract_grid::<SingleElementGrid<f32, CiarletElement<f32>>>(grid))
                    .topology_dim(),
                DType::F64 => (*extract_grid::<SingleElementGrid<f64, CiarletElement<f64>>>(grid))
                    .topology_dim(),
            },
        }
    }

    #[no_mangle]
    pub unsafe extern "C" fn grid_gdim(grid: *const GridWrapper) -> usize {
        match (*grid).gtype {
            GridType::SerialSingleElementGrid => match (*grid).dtype {
                DType::F32 => (*extract_grid::<SingleElementGrid<f32, CiarletElement<f32>>>(grid))
                    .geometry_dim(),
                DType::F64 => (*extract_grid::<SingleElementGrid<f64, CiarletElement<f64>>>(grid))
                    .geometry_dim(),
            },
        }
    }

    #[no_mangle]
    pub unsafe extern "C" fn grid_entity_count(grid: *const GridWrapper, entity_type: u8) -> usize {
        match (*grid).gtype {
            GridType::SerialSingleElementGrid => match (*grid).dtype {
                DType::F32 => (*extract_grid::<SingleElementGrid<f32, CiarletElement<f32>>>(grid))
                    .entity_count(ReferenceCellType::from(entity_type).unwrap()),
                DType::F64 => (*extract_grid::<SingleElementGrid<f64, CiarletElement<f64>>>(grid))
                    .entity_count(ReferenceCellType::from(entity_type).unwrap()),
            },
        }
    }

    unsafe fn grid_entity_internal<T: Grid>(
        grid: *const GridWrapper,
        dim: usize,
        local_index: usize,
        etype: EntityType,
    ) -> *const EntityWrapper {
        let entity = EntityWrapper {
            entity: Box::into_raw(Box::new(
                (*extract_grid::<T>(grid)).entity(dim, local_index).unwrap(),
            )) as *const c_void,
            dtype: (*grid).dtype,
            etype,
        };
        Box::into_raw(Box::new(entity))
    }

    #[no_mangle]
    pub unsafe extern "C" fn grid_entity(
        grid: *const GridWrapper,
        dim: usize,
        local_index: usize,
    ) -> *const EntityWrapper {
        match (*grid).gtype {
            GridType::SerialSingleElementGrid => match (*grid).dtype {
                DType::F32 => grid_entity_internal::<SingleElementGrid<f32, CiarletElement<f32>>>(
                    grid,
                    dim,
                    local_index,
                    EntityType::SingleElementGridEntity,
                ),
                DType::F64 => grid_entity_internal::<SingleElementGrid<f64, CiarletElement<f64>>>(
                    grid,
                    dim,
                    local_index,
                    EntityType::SingleElementGridEntity,
                ),
            },
        }
    }

    unsafe fn grid_entity_from_id_internal<T: Grid>(
        grid: *const GridWrapper,
        dim: usize,
        id: usize,
        etype: EntityType,
    ) -> *const EntityWrapper {
        let entity = EntityWrapper {
            entity: Box::into_raw(Box::new(
                (*extract_grid::<T>(grid)).entity_from_id(dim, id).unwrap(),
            )) as *const c_void,
            dtype: (*grid).dtype,
            etype,
        };
        Box::into_raw(Box::new(entity))
    }

    #[no_mangle]
    pub unsafe extern "C" fn grid_entity_from_id(
        grid: *const GridWrapper,
        dim: usize,
        id: usize,
    ) -> *const EntityWrapper {
        match (*grid).gtype {
            GridType::SerialSingleElementGrid => match (*grid).dtype {
                DType::F32 => grid_entity_from_id_internal::<
                    SingleElementGrid<f32, CiarletElement<f32>>,
                >(grid, dim, id, EntityType::SingleElementGridEntity),
                DType::F64 => grid_entity_from_id_internal::<
                    SingleElementGrid<f64, CiarletElement<f64>>,
                >(grid, dim, id, EntityType::SingleElementGridEntity),
            },
        }
    }

    #[no_mangle]
    pub unsafe extern "C" fn grid_entity_types_size(grid: *const GridWrapper, dim: usize) -> usize {
        match (*grid).gtype {
            GridType::SerialSingleElementGrid => match (*grid).dtype {
                DType::F32 => (*extract_grid::<SingleElementGrid<f32, CiarletElement<f32>>>(grid))
                    .entity_types(dim),
                DType::F64 => (*extract_grid::<SingleElementGrid<f64, CiarletElement<f64>>>(grid))
                    .entity_types(dim),
            },
        }
        .len()
    }

    #[no_mangle]
    pub unsafe extern "C" fn grid_entity_types(
        grid: *const GridWrapper,
        dim: usize,
        types: *mut u8,
    ) {
        for (i, t) in match (*grid).gtype {
            GridType::SerialSingleElementGrid => match (*grid).dtype {
                DType::F32 => (*extract_grid::<SingleElementGrid<f32, CiarletElement<f32>>>(grid))
                    .entity_types(dim),
                DType::F64 => (*extract_grid::<SingleElementGrid<f64, CiarletElement<f64>>>(grid))
                    .entity_types(dim),
            },
        }
        .iter()
        .enumerate()
        {
            *types.add(i) = *t as u8;
        }
    }

    #[no_mangle]
    pub unsafe extern "C" fn grid_dtype(grid: *const GridWrapper) -> u8 {
        (*grid).dtype as u8
    }
}

mod entity {
    use super::geometry::{GeometryType, GeometryWrapper};
    use super::topology::{TopologyType, TopologyWrapper};
    use super::DType;
    use crate::{grid::serial::SingleElementGridEntity, traits::Entity, types::Ownership};
    use ndelement::ciarlet::CiarletElement;
    use std::ffi::c_void;

    #[derive(Debug, PartialEq, Clone, Copy)]
    #[repr(u8)]
    pub enum EntityType {
        SingleElementGridEntity = 0,
    }

    #[repr(C)]
    pub struct EntityWrapper {
        pub entity: *const c_void,
        pub dtype: DType,
        pub etype: EntityType,
    }

    impl Drop for EntityWrapper {
        fn drop(&mut self) {
            let Self {
                entity,
                dtype,
                etype,
            } = self;
            match etype {
                EntityType::SingleElementGridEntity => match dtype {
                    DType::F32 => drop(unsafe {
                        Box::from_raw(
                            *entity as *mut SingleElementGridEntity<f32, CiarletElement<f32>>,
                        )
                    }),
                    DType::F64 => drop(unsafe {
                        Box::from_raw(
                            *entity as *mut SingleElementGridEntity<f64, CiarletElement<f64>>,
                        )
                    }),
                },
            }
        }
    }

    #[no_mangle]
    pub unsafe extern "C" fn free_entity(e: *mut EntityWrapper) {
        assert!(!e.is_null());
        unsafe { drop(Box::from_raw(e)) }
    }

    unsafe fn extract_entity<E: Entity>(entity: *const EntityWrapper) -> *const E {
        (*entity).entity as *const E
    }

    #[no_mangle]
    pub unsafe extern "C" fn entity_local_index(entity: *const EntityWrapper) -> usize {
        match (*entity).etype {
            EntityType::SingleElementGridEntity => match (*entity).dtype {
                DType::F32 => {
                    (*extract_entity::<SingleElementGridEntity<f32, CiarletElement<f32>>>(entity))
                        .local_index()
                }
                DType::F64 => {
                    (*extract_entity::<SingleElementGridEntity<f64, CiarletElement<f64>>>(entity))
                        .local_index()
                }
            },
        }
    }

    #[no_mangle]
    pub unsafe extern "C" fn entity_global_index(entity: *const EntityWrapper) -> usize {
        match (*entity).etype {
            EntityType::SingleElementGridEntity => match (*entity).dtype {
                DType::F32 => {
                    (*extract_entity::<SingleElementGridEntity<f32, CiarletElement<f32>>>(entity))
                        .global_index()
                }
                DType::F64 => {
                    (*extract_entity::<SingleElementGridEntity<f64, CiarletElement<f64>>>(entity))
                        .global_index()
                }
            },
        }
    }

    #[no_mangle]
    pub unsafe extern "C" fn entity_entity_type(entity: *const EntityWrapper) -> u8 {
        match (*entity).etype {
            EntityType::SingleElementGridEntity => match (*entity).dtype {
                DType::F32 => {
                    (*extract_entity::<SingleElementGridEntity<f32, CiarletElement<f32>>>(entity))
                        .entity_type() as u8
                }

                DType::F64 => {
                    (*extract_entity::<SingleElementGridEntity<f64, CiarletElement<f64>>>(entity))
                        .entity_type() as u8
                }
            },
        }
    }
    #[no_mangle]
    pub unsafe extern "C" fn entity_has_id(entity: *const EntityWrapper) -> bool {
        match (*entity).etype {
            EntityType::SingleElementGridEntity => match (*entity).dtype {
                DType::F32 => {
                    (*extract_entity::<SingleElementGridEntity<f32, CiarletElement<f32>>>(entity))
                        .id()
                }
                DType::F64 => {
                    (*extract_entity::<SingleElementGridEntity<f64, CiarletElement<f64>>>(entity))
                        .id()
                }
            },
        }
        .is_some()
    }

    #[no_mangle]
    pub unsafe extern "C" fn entity_id(entity: *const EntityWrapper) -> usize {
        match (*entity).etype {
            EntityType::SingleElementGridEntity => match (*entity).dtype {
                DType::F32 => {
                    (*extract_entity::<SingleElementGridEntity<f32, CiarletElement<f32>>>(entity))
                        .id()
                }
                DType::F64 => {
                    (*extract_entity::<SingleElementGridEntity<f64, CiarletElement<f64>>>(entity))
                        .id()
                }
            },
        }
        .unwrap()
    }

    #[no_mangle]
    pub unsafe extern "C" fn entity_is_owned(entity: *const EntityWrapper) -> bool {
        let ownership = match (*entity).etype {
            EntityType::SingleElementGridEntity => match (*entity).dtype {
                DType::F32 => {
                    (*extract_entity::<SingleElementGridEntity<f32, CiarletElement<f32>>>(entity))
                        .ownership()
                }
                DType::F64 => {
                    (*extract_entity::<SingleElementGridEntity<f64, CiarletElement<f64>>>(entity))
                        .ownership()
                }
            },
        };

        ownership == Ownership::Owned
    }

    #[no_mangle]
    pub unsafe extern "C" fn entity_ownership_process(entity: *const EntityWrapper) -> usize {
        let ownership = match (*entity).etype {
            EntityType::SingleElementGridEntity => match (*entity).dtype {
                DType::F32 => {
                    (*extract_entity::<SingleElementGridEntity<f32, CiarletElement<f32>>>(entity))
                        .ownership()
                }
                DType::F64 => {
                    (*extract_entity::<SingleElementGridEntity<f64, CiarletElement<f64>>>(entity))
                        .ownership()
                }
            },
        };

        if let Ownership::Ghost(process, _index) = ownership {
            process
        } else {
            panic!();
        }
    }

    #[no_mangle]
    pub unsafe extern "C" fn entity_ownership_index(entity: *const EntityWrapper) -> usize {
        let ownership = match (*entity).etype {
            EntityType::SingleElementGridEntity => match (*entity).dtype {
                DType::F32 => {
                    (*extract_entity::<SingleElementGridEntity<f32, CiarletElement<f32>>>(entity))
                        .ownership()
                }
                DType::F64 => {
                    (*extract_entity::<SingleElementGridEntity<f64, CiarletElement<f64>>>(entity))
                        .ownership()
                }
            },
        };

        if let Ownership::Ghost(_process, index) = ownership {
            index
        } else {
            panic!();
        }
    }

    unsafe fn entity_topology_internal<E: Entity>(
        entity: *const EntityWrapper,
        ttype: TopologyType,
    ) -> *const TopologyWrapper {
        let topology = TopologyWrapper {
            topology: Box::into_raw(Box::new((*extract_entity::<E>(entity)).topology()))
                as *const c_void,
            ttype,
        };
        Box::into_raw(Box::new(topology))
    }

    #[no_mangle]
    pub unsafe extern "C" fn entity_topology(
        entity: *const EntityWrapper,
    ) -> *const TopologyWrapper {
        match (*entity).etype {
            EntityType::SingleElementGridEntity => match (*entity).dtype {
                DType::F32 => entity_topology_internal::<
                    SingleElementGridEntity<f32, CiarletElement<f32>>,
                >(entity, TopologyType::SingleTypeEntityTopology),
                DType::F64 => entity_topology_internal::<
                    SingleElementGridEntity<f64, CiarletElement<f64>>,
                >(entity, TopologyType::SingleTypeEntityTopology),
            },
        }
    }

    unsafe fn entity_geometry_internal<E: Entity>(
        entity: *const EntityWrapper,
        gtype: GeometryType,
        dtype: DType,
    ) -> *const GeometryWrapper {
        let geometry = GeometryWrapper {
            geometry: Box::into_raw(Box::new((*extract_entity::<E>(entity)).geometry()))
                as *const c_void,
            gtype,
            dtype,
        };
        Box::into_raw(Box::new(geometry))
    }

    #[no_mangle]
    pub unsafe extern "C" fn entity_geometry(
        entity: *const EntityWrapper,
    ) -> *const GeometryWrapper {
        match (*entity).etype {
            EntityType::SingleElementGridEntity => match (*entity).dtype {
                DType::F32 => {
                    entity_geometry_internal::<SingleElementGridEntity<f32, CiarletElement<f32>>>(
                        entity,
                        GeometryType::SingleElementEntityGeometry,
                        (*entity).dtype,
                    )
                }
                DType::F64 => {
                    entity_geometry_internal::<SingleElementGridEntity<f64, CiarletElement<f64>>>(
                        entity,
                        GeometryType::SingleElementEntityGeometry,
                        (*entity).dtype,
                    )
                }
            },
        }
    }
    #[no_mangle]
    pub unsafe extern "C" fn entity_dtype(entity: *const EntityWrapper) -> u8 {
        (*entity).dtype as u8
    }
}

mod topology {
    use crate::{topology::serial::SingleTypeEntityTopology, traits::Topology};
    use std::ffi::c_void;

    #[derive(Debug, PartialEq, Clone, Copy)]
    #[repr(u8)]
    pub enum TopologyType {
        SingleTypeEntityTopology = 0,
    }

    #[repr(C)]
    pub struct TopologyWrapper {
        pub topology: *const c_void,
        pub ttype: TopologyType,
    }

    impl Drop for TopologyWrapper {
        fn drop(&mut self) {
            let Self { topology, ttype } = self;
            match ttype {
                TopologyType::SingleTypeEntityTopology => {
                    drop(unsafe { Box::from_raw(*topology as *mut SingleTypeEntityTopology) })
                }
            }
        }
    }

    #[no_mangle]
    pub unsafe extern "C" fn free_topology(t: *mut TopologyWrapper) {
        assert!(!t.is_null());
        unsafe { drop(Box::from_raw(t)) }
    }

    unsafe fn extract_topology<T: Topology>(topology: *const TopologyWrapper) -> *const T {
        (*topology).topology as *const T
    }

    #[no_mangle]
    pub unsafe extern "C" fn topology_sub_entity(
        topology: *const TopologyWrapper,
        dim: usize,
        index: usize,
    ) -> usize {
        match (*topology).ttype {
            TopologyType::SingleTypeEntityTopology => {
                (*extract_topology::<SingleTypeEntityTopology>(topology)).sub_entity(dim, index)
            }
        }
    }
}

mod geometry {
    use super::DType;
    use crate::geometry::SingleElementEntityGeometry;
    use ndelement::ciarlet::CiarletElement;
    use std::ffi::c_void;

    #[derive(Debug, PartialEq, Clone, Copy)]
    #[repr(u8)]
    pub enum GeometryType {
        SingleElementEntityGeometry = 0,
    }

    #[repr(C)]
    pub struct GeometryWrapper {
        pub geometry: *const c_void,
        pub dtype: DType,
        pub gtype: GeometryType,
    }

    impl Drop for GeometryWrapper {
        fn drop(&mut self) {
            let Self {
                geometry,
                dtype,
                gtype,
            } = self;
            match gtype {
                GeometryType::SingleElementEntityGeometry => match dtype {
                    DType::F32 => drop(unsafe {
                        Box::from_raw(
                            *geometry as *mut SingleElementEntityGeometry<f32, CiarletElement<f32>>,
                        )
                    }),
                    DType::F64 => drop(unsafe {
                        Box::from_raw(
                            *geometry as *mut SingleElementEntityGeometry<f64, CiarletElement<f64>>,
                        )
                    }),
                },
            }
        }
    }

    #[no_mangle]
    pub unsafe extern "C" fn free_geometry(g: *mut GeometryWrapper) {
        assert!(!g.is_null());
        unsafe { drop(Box::from_raw(g)) }
    }

    //unsafe fn extract_geometry<G: Geometry>(geometry: *const GeometryWrapper) -> *const G {
    //    (*geometry).geometry as *const G
    //}
}

mod shapes {
    use super::grid::{GridType, GridWrapper};
    use super::DType;
    use crate::shapes::regular_sphere;
    use crate::types::RealScalar;
    use std::ffi::c_void;

    unsafe fn shapes_regular_sphere<T: RealScalar>(level: u32, dtype: DType) -> *const GridWrapper {
        let grid = Box::into_raw(Box::new(regular_sphere::<T>(level))) as *const c_void;
        Box::into_raw(Box::new(GridWrapper {
            grid,
            dtype,
            gtype: GridType::SerialSingleElementGrid,
        }))
    }
    #[no_mangle]
    pub unsafe extern "C" fn shapes_regular_sphere_f64(level: u32) -> *const GridWrapper {
        shapes_regular_sphere::<f64>(level, DType::F64)
    }
    #[no_mangle]
    pub unsafe extern "C" fn shapes_regular_sphere_f32(level: u32) -> *const GridWrapper {
        shapes_regular_sphere::<f32>(level, DType::F32)
    }
}
