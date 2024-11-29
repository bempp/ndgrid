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
    use crate::{traits::Grid, types::RealScalar, SingleElementGrid, SingleElementGridBorrowed};
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
        pub fn single_element_grid_impl<T: RealScalar>(
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
            DType::F32 => single_element_grid_impl(
                coordinates as *const f32,
                npts,
                gdim,
                cells,
                ncells,
                cell_type,
                geometry_degree,
            ),
            DType::F64 => single_element_grid_impl(
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

    #[repr(u8)]
    #[derive(Debug)]
    pub enum GridType {
        SingleElementGrid,
        SingleElementGridBorrowed,
    }

    #[no_mangle]
    pub unsafe extern "C" fn grid_type(grid: *mut GridT) -> GridType {
        #[concretise_types(
            gen_type(name = "dtype", replace_with = ["f32", "f64"]),
            field(arg = 1, name = "wrap", wrapper = "GridT", replace_with = ["SingleElementGrid<{{dtype}}, CiarletElement<{{dtype}}>>", "SingleElementGridBorrowed<{{dtype}}, CiarletElement<{{dtype}}>>"]),
        )]
        pub unsafe fn grid_type_internal<T: RealScalar, G: Grid<T = T>>(
            grid_t: *mut GridT,
            _grid: &G,
        ) -> GridType {
            if grid_t
                .as_ref()
                .unwrap()
                .inner()
                .downcast_ref::<SingleElementGrid<T, CiarletElement<T>>>()
                .is_some()
            {
                GridType::SingleElementGrid
            } else if grid_t
                .as_ref()
                .unwrap()
                .inner()
                .downcast_ref::<SingleElementGridBorrowed<T, CiarletElement<T>>>()
                .is_some()
            {
                GridType::SingleElementGridBorrowed
            } else {
                {
                    panic!("Unknown type.");
                };
            }
        }
        grid_type_internal(grid, grid)
    }

    pub mod single_element_grid {
        use super::super::{grid_t_create, grid_t_unwrap, GridT};
        use crate::{
            traits::Grid, types::RealScalar, SingleElementGrid, SingleElementGridBorrowed,
        };
        use c_api_tools::{concretise_types, DType, DTypeIdentifier};
        use core::ffi::c_void;
        use ndelement::{
            ciarlet::{CiarletElement, LagrangeElementFamily},
            traits::ElementFamily,
            types::{Continuity, ReferenceCellType},
        };
        use rlst::{rlst_array_from_slice2, RawAccess, Shape};
        use std::slice::from_raw_parts;

        pub struct InternalDataContainer {
            _ptr: Box<dyn std::any::Any>,
        }

        /// Free the instance of the wrapper.
        #[no_mangle]
        pub unsafe extern "C" fn internal_data_container_free(ptr: *mut InternalDataContainer) {
            if ptr.is_null() {
                return;
            }
            unsafe {
                drop(Box::from_raw(ptr));
            }
        }

        #[repr(C)]
        pub struct SingleElementGridCData {
            internal_storage: *const InternalDataContainer,
            tdim: usize,
            id_sizes: *const usize,
            id_pointers: *const *const usize,
            entity_types: *const ReferenceCellType,
            entity_counts: *const usize,
            downward_connectivity: *const *const *const usize,
            downward_connectivity_shape0: *const *const usize,
            upward_connectivity: *const *const *const *const usize,
            upward_connectivity_lens: *const *const *const usize,
            points: *const c_void,
            gdim: usize,
            npoints: usize,
            dtype: DType,
            cells: *const usize,
            points_per_cell: usize,
            ncells: usize,
            geometry_degree: usize,
        }

        #[allow(clippy::type_complexity)]
        struct SingleElementGridInternalData {
            _id_sizes: Vec<usize>,
            _id_pointers: Vec<*const usize>,
            _dc: (Vec<Vec<*const usize>>, Vec<*const *const usize>),
            _dcsh: (Vec<Vec<usize>>, Vec<*const usize>),
            _uc: (
                Vec<Vec<Vec<*const usize>>>,
                Vec<Vec<*const *const usize>>,
                Vec<*const *const *const usize>,
            ),
            _ucl: (
                Vec<Vec<Vec<usize>>>,
                Vec<Vec<*const usize>>,
                Vec<*const *const usize>,
            ),
        }

        #[concretise_types(
            gen_type(name = "dtype", replace_with = ["f32", "f64"]),
            field(arg = 0, name = "wrap", wrapper = "GridT", replace_with = ["SingleElementGrid<{{dtype}}, CiarletElement<{{dtype}}>>"]),
        )]
        pub fn single_element_grid_cdata<T: RealScalar + DTypeIdentifier>(
            grid: &SingleElementGrid<T, CiarletElement<T>>,
        ) -> SingleElementGridCData {
            let tdim = grid.topology_dim();
            let gdim = grid.geometry_dim();
            let dtype = <T as DTypeIdentifier>::dtype();

            let id_sizes_data = grid
                .internal_topology()
                .ids
                .iter()
                .map(|i| if let Some(j) = i { j.len() } else { 0 })
                .collect::<Vec<_>>();
            let id_sizes = id_sizes_data.as_ptr();

            let null_data = [];
            let id_pointers_data = grid
                .internal_topology()
                .ids
                .iter()
                .map(|i| {
                    if let Some(j) = i {
                        j.as_ptr()
                    } else {
                        null_data.as_ptr()
                    }
                })
                .collect::<Vec<_>>();
            let id_pointers = id_pointers_data.as_ptr();

            let entity_types = grid.internal_topology().entity_types().as_ptr();
            let entity_counts = grid.internal_topology().entity_counts().as_ptr();

            let dc0 = grid
                .internal_topology()
                .downward_connectivity()
                .iter()
                .map(|i| i.iter().map(|j| j.data().as_ptr()).collect::<Vec<_>>())
                .collect::<Vec<_>>();
            let dc1 = dc0.iter().map(|i| i.as_ptr()).collect::<Vec<_>>();
            let dc = (dc0, dc1);
            let downward_connectivity = dc.1.as_ptr();

            let dcsh0 = grid
                .internal_topology()
                .downward_connectivity()
                .iter()
                .map(|i| i.iter().map(|j| j.shape()[0]).collect::<Vec<_>>())
                .collect::<Vec<_>>();
            let dcsh1 = dcsh0.iter().map(|i| i.as_ptr()).collect::<Vec<_>>();
            let dcsh = (dcsh0, dcsh1);
            let downward_connectivity_shape0 = dcsh.1.as_ptr();

            let uc0 = grid
                .internal_topology()
                .upward_connectivity()
                .iter()
                .map(|i| {
                    i.iter()
                        .map(|j| j.iter().map(|k| k.as_ptr()).collect::<Vec<_>>())
                        .collect::<Vec<_>>()
                })
                .collect::<Vec<_>>();
            let uc1 = uc0
                .iter()
                .map(|i| i.iter().map(|j| j.as_ptr()).collect::<Vec<_>>())
                .collect::<Vec<_>>();
            let uc2 = uc1.iter().map(|i| i.as_ptr()).collect::<Vec<_>>();
            let uc = (uc0, uc1, uc2);
            let upward_connectivity = uc.2.as_ptr();
            let ucl0 = grid
                .internal_topology()
                .upward_connectivity()
                .iter()
                .map(|i| {
                    i.iter()
                        .map(|j| j.iter().map(|k| k.len()).collect::<Vec<_>>())
                        .collect::<Vec<_>>()
                })
                .collect::<Vec<_>>();
            let ucl1 = ucl0
                .iter()
                .map(|i| i.iter().map(|j| j.as_ptr()).collect::<Vec<_>>())
                .collect::<Vec<_>>();
            let ucl2 = ucl1.iter().map(|i| i.as_ptr()).collect::<Vec<_>>();
            let ucl = (ucl0, ucl1, ucl2);
            let upward_connectivity_lens = ucl.2.as_ptr();

            let points = grid.internal_geometry().points().data().as_ptr() as *const c_void;
            let npoints = grid.internal_geometry().points().shape()[1];
            let cells = grid.internal_geometry().cells().data().as_ptr();
            let [points_per_cell, ncells] = grid.internal_geometry().cells().shape();
            let geometry_degree = grid.internal_geometry().element().degree();

            let obj = InternalDataContainer { _ptr: Box::new(()) };
            let internal_storage = Box::into_raw(Box::new(obj));
            unsafe {
                (*internal_storage)._ptr = Box::new(SingleElementGridInternalData {
                    _id_sizes: id_sizes_data,
                    _id_pointers: id_pointers_data,
                    _dc: dc,
                    _dcsh: dcsh,
                    _uc: uc,
                    _ucl: ucl,
                });
            }
            SingleElementGridCData {
                internal_storage,
                tdim,
                id_sizes,
                id_pointers,
                entity_types,
                entity_counts,
                downward_connectivity,
                downward_connectivity_shape0,
                upward_connectivity,
                upward_connectivity_lens,
                points,
                gdim,
                npoints,
                dtype,
                cells,
                points_per_cell,
                ncells,
                geometry_degree,
            }
        }

        #[no_mangle]
        pub unsafe extern "C" fn single_element_grid_borrowed_create(
            tdim: usize,
            id_sizes: *const usize,
            id_pointers: *const *const usize,
            entity_types: *const ReferenceCellType,
            entity_counts: *const usize,
            downward_connectivity: *const *const *const usize,
            downward_connectivity_shape0: *const *const usize,
            upward_connectivity: *const *const *const *const usize,
            upward_connectivity_lens: *const *const *const usize,
            points: *const c_void,
            gdim: usize,
            npoints: usize,
            dtype: DType,
            cells: *const usize,
            points_per_cell: usize,
            ncells: usize,
            geometry_degree: usize,
        ) -> *mut GridT {
            let wrapper = grid_t_create();
            let inner = unsafe { grid_t_unwrap(wrapper) }.unwrap();

            let ids = (0..=tdim)
                .map(|d| {
                    if *id_sizes.add(d) == 0 {
                        None
                    } else {
                        Some(from_raw_parts(*id_pointers.add(d), *id_sizes.add(d)))
                    }
                })
                .collect::<Vec<_>>();
            let entity_types = from_raw_parts(entity_types, tdim + 1);
            let entity_counts = from_raw_parts(entity_counts, tdim + 1);
            let downward_connectivity = (0..tdim + 1)
                .map(|d| {
                    (0..d + 1)
                        .map(|i| {
                            let shape = [
                                *(*downward_connectivity_shape0.add(d)).add(i),
                                entity_counts[i],
                            ];
                            rlst_array_from_slice2!(
                                from_raw_parts(
                                    *(*downward_connectivity.add(d)).add(i),
                                    shape[0] * shape[1]
                                ),
                                shape
                            )
                        })
                        .collect::<Vec<_>>()
                })
                .collect::<Vec<_>>();

            let upward_connectivity = (0..tdim)
                .map(|d| {
                    (0..tdim - d)
                        .map(|i| {
                            (0..entity_counts[d])
                                .map(|j| {
                                    from_raw_parts(
                                        *(*(*upward_connectivity.add(d)).add(i)).add(j),
                                        *(*(*upward_connectivity_lens.add(d)).add(i)).add(j),
                                    )
                                })
                                .collect::<Vec<_>>()
                        })
                        .collect::<Vec<_>>()
                })
                .collect::<Vec<_>>();
            let cells = rlst_array_from_slice2!(
                from_raw_parts(cells, points_per_cell * ncells),
                [points_per_cell, ncells]
            );

            match dtype {
                DType::F32 => {
                    let points = rlst_array_from_slice2!(
                        from_raw_parts(points as *const f32, gdim * npoints),
                        [gdim, npoints]
                    );
                    let family = LagrangeElementFamily::new(geometry_degree, Continuity::Standard);
                    let elements = entity_types
                        .iter()
                        .skip(1)
                        .map(|t| family.element(*t))
                        .collect::<Vec<_>>();
                    *inner = Box::new(SingleElementGridBorrowed::new(
                        tdim,
                        ids,
                        entity_types,
                        entity_counts,
                        downward_connectivity,
                        upward_connectivity,
                        points,
                        cells,
                        elements,
                    ));
                }
                DType::F64 => {
                    let points = rlst_array_from_slice2!(
                        from_raw_parts(points as *const f64, gdim * npoints),
                        [gdim, npoints]
                    );
                    let family = LagrangeElementFamily::new(geometry_degree, Continuity::Standard);
                    let elements = entity_types
                        .iter()
                        .skip(1)
                        .map(|t| family.element(*t))
                        .collect::<Vec<_>>();
                    *inner = Box::new(SingleElementGridBorrowed::new(
                        tdim,
                        ids,
                        entity_types,
                        entity_counts,
                        downward_connectivity,
                        upward_connectivity,
                        points,
                        cells,
                        elements,
                    ));
                }
                _ => {
                    panic!("Unsupported dtype");
                }
            }

            wrapper
        }
    }

    #[concretise_types(
        gen_type(name = "dtype", replace_with = ["f32", "f64"]),
        field(arg = 0, name = "wrap", wrapper = "GridT", replace_with = ["SingleElementGrid<{{dtype}}, CiarletElement<{{dtype}}>>", "SingleElementGridBorrowed<{{dtype}}, CiarletElement<{{dtype}}>>"]),
    )]
    pub fn grid_tdim<G: Grid>(grid: &G) -> usize {
        grid.topology_dim()
    }

    #[concretise_types(
        gen_type(name = "dtype", replace_with = ["f32", "f64"]),
        field(arg = 0, name = "wrap", wrapper = "GridT", replace_with = ["SingleElementGrid<{{dtype}}, CiarletElement<{{dtype}}>>", "SingleElementGridBorrowed<{{dtype}}, CiarletElement<{{dtype}}>>"]),
    )]
    pub fn grid_gdim<G: Grid>(grid: &G) -> usize {
        grid.geometry_dim()
    }

    #[concretise_types(
        gen_type(name = "dtype", replace_with = ["f32", "f64"]),
        field(arg = 0, name = "wrap", wrapper = "GridT", replace_with = ["SingleElementGrid<{{dtype}}, CiarletElement<{{dtype}}>>", "SingleElementGridBorrowed<{{dtype}}, CiarletElement<{{dtype}}>>"]),
    )]
    pub fn grid_entity_count<G: Grid<EntityDescriptor = ReferenceCellType>>(
        grid: &G,
        entity_type: ReferenceCellType,
    ) -> usize {
        grid.entity_count(entity_type)
    }

    #[concretise_types(
        gen_type(name = "dtype", replace_with = ["f32", "f64"]),
        field(arg = 0, name = "wrap", wrapper = "GridT", replace_with = ["SingleElementGrid<{{dtype}}, CiarletElement<{{dtype}}>>", "SingleElementGridBorrowed<{{dtype}}, CiarletElement<{{dtype}}>>"]),
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
        field(arg = 0, name = "wrap", wrapper = "GridT", replace_with = ["SingleElementGrid<{{dtype}}, CiarletElement<{{dtype}}>>", "SingleElementGridBorrowed<{{dtype}}, CiarletElement<{{dtype}}>>"]),
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
        field(arg = 0, name = "wrap", wrapper = "GridT", replace_with = ["SingleElementGrid<{{dtype}}, CiarletElement<{{dtype}}>>", "SingleElementGridBorrowed<{{dtype}}, CiarletElement<{{dtype}}>>"]),
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
        field(arg = 0, name = "wrap", wrapper = "GridT", replace_with = ["SingleElementGrid<{{dtype}}, CiarletElement<{{dtype}}>>", "SingleElementGridBorrowed<{{dtype}}, CiarletElement<{{dtype}}>>"]),
    )]
    pub fn grid_entity_types_size<G: Grid>(grid: &G, dim: usize) -> usize {
        grid.entity_types(dim).len()
    }

    #[concretise_types(
        gen_type(name = "dtype", replace_with = ["f32", "f64"]),
        field(arg = 0, name = "wrap", wrapper = "GridT", replace_with = ["SingleElementGrid<{{dtype}}, CiarletElement<{{dtype}}>>", "SingleElementGridBorrowed<{{dtype}}, CiarletElement<{{dtype}}>>"]),
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
        field(arg = 0, name = "wrap", wrapper = "GridT", replace_with = ["SingleElementGrid<{{dtype}}, CiarletElement<{{dtype}}>>", "SingleElementGridBorrowed<{{dtype}}, CiarletElement<{{dtype}}>>"]),
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
        grid::serial::{SingleElementGridEntity, SingleElementGridEntityBorrowed},
        traits::Entity,
        types::{Ownership, RealScalar},
    };
    use c_api_tools::{concretise_types, DType, DTypeIdentifier};
    use ndelement::{ciarlet::CiarletElement, types::ReferenceCellType};

    #[concretise_types(
        gen_type(name = "dtype", replace_with = ["f32", "f64"]),
        field(arg = 0, name = "wrap", wrapper = "EntityT", replace_with = ["SingleElementGridEntity<{{dtype}}, CiarletElement<{{dtype}}>>", "SingleElementGridEntityBorrowed<{{dtype}}, CiarletElement<{{dtype}}>>"]),
    )]
    pub fn entity_local_index<E: Entity>(entity: &E) -> usize {
        entity.local_index()
    }

    #[concretise_types(
        gen_type(name = "dtype", replace_with = ["f32", "f64"]),
        field(arg = 0, name = "wrap", wrapper = "EntityT", replace_with = ["SingleElementGridEntity<{{dtype}}, CiarletElement<{{dtype}}>>", "SingleElementGridEntityBorrowed<{{dtype}}, CiarletElement<{{dtype}}>>"]),
    )]
    pub fn entity_global_index<E: Entity>(entity: &E) -> usize {
        entity.global_index()
    }

    #[concretise_types(
        gen_type(name = "dtype", replace_with = ["f32", "f64"]),
        field(arg = 0, name = "wrap", wrapper = "EntityT", replace_with = ["SingleElementGridEntity<{{dtype}}, CiarletElement<{{dtype}}>>", "SingleElementGridEntityBorrowed<{{dtype}}, CiarletElement<{{dtype}}>>"]),
    )]
    pub fn entity_entity_type<E: Entity<EntityDescriptor = ReferenceCellType>>(
        entity: &E,
    ) -> ReferenceCellType {
        entity.entity_type()
    }

    #[concretise_types(
        gen_type(name = "dtype", replace_with = ["f32", "f64"]),
        field(arg = 0, name = "wrap", wrapper = "EntityT", replace_with = ["SingleElementGridEntity<{{dtype}}, CiarletElement<{{dtype}}>>", "SingleElementGridEntityBorrowed<{{dtype}}, CiarletElement<{{dtype}}>>"]),
    )]
    pub fn entity_has_id<E: Entity>(entity: &E) -> bool {
        entity.id().is_some()
    }

    #[concretise_types(
        gen_type(name = "dtype", replace_with = ["f32", "f64"]),
        field(arg = 0, name = "wrap", wrapper = "EntityT", replace_with = ["SingleElementGridEntity<{{dtype}}, CiarletElement<{{dtype}}>>", "SingleElementGridEntityBorrowed<{{dtype}}, CiarletElement<{{dtype}}>>"]),
    )]
    pub fn entity_id<E: Entity>(entity: &E) -> usize {
        entity.id().unwrap()
    }

    #[concretise_types(
        gen_type(name = "dtype", replace_with = ["f32", "f64"]),
        field(arg = 0, name = "wrap", wrapper = "EntityT", replace_with = ["SingleElementGridEntity<{{dtype}}, CiarletElement<{{dtype}}>>", "SingleElementGridEntityBorrowed<{{dtype}}, CiarletElement<{{dtype}}>>"]),
    )]
    pub fn entity_is_owned<E: Entity>(entity: &E) -> bool {
        entity.ownership() == Ownership::Owned
    }

    #[concretise_types(
        gen_type(name = "dtype", replace_with = ["f32", "f64"]),
        field(arg = 0, name = "wrap", wrapper = "EntityT", replace_with = ["SingleElementGridEntity<{{dtype}}, CiarletElement<{{dtype}}>>", "SingleElementGridEntityBorrowed<{{dtype}}, CiarletElement<{{dtype}}>>"]),
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
        field(arg = 0, name = "wrap", wrapper = "EntityT", replace_with = ["SingleElementGridEntity<{{dtype}}, CiarletElement<{{dtype}}>>", "SingleElementGridEntityBorrowed<{{dtype}}, CiarletElement<{{dtype}}>>"]),
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
        field(arg = 0, name = "wrap", wrapper = "EntityT", replace_with = ["SingleElementGridEntity<{{dtype}}, CiarletElement<{{dtype}}>>", "SingleElementGridEntityBorrowed<{{dtype}}, CiarletElement<{{dtype}}>>"]),
    )]
    pub fn entity_topology<E: Entity>(entity: &'static E) -> *mut TopologyT {
        let wrapper = topology_t_create();
        let inner = unsafe { topology_t_unwrap(wrapper) }.unwrap();
        *inner = Box::new(entity.topology());
        wrapper
    }

    #[concretise_types(
        gen_type(name = "dtype", replace_with = ["f32", "f64"]),
        field(arg = 0, name = "wrap", wrapper = "EntityT", replace_with = ["SingleElementGridEntity<{{dtype}}, CiarletElement<{{dtype}}>>", "SingleElementGridEntityBorrowed<{{dtype}}, CiarletElement<{{dtype}}>>"]),
    )]
    pub fn entity_geometry<E: Entity>(entity: &'static E) -> *mut GeometryT {
        let wrapper = geometry_t_create();
        let inner = unsafe { geometry_t_unwrap(wrapper) }.unwrap();
        *inner = Box::new(entity.geometry());
        wrapper
    }

    #[concretise_types(
        gen_type(name = "dtype", replace_with = ["f32", "f64"]),
        field(arg = 0, name = "wrap", wrapper = "EntityT", replace_with = ["SingleElementGridEntity<{{dtype}}, CiarletElement<{{dtype}}>>", "SingleElementGridEntityBorrowed<{{dtype}}, CiarletElement<{{dtype}}>>"]),
    )]
    pub fn entity_dtype<T: RealScalar + DTypeIdentifier, E: Entity<T = T>>(_entity: &E) -> DType {
        <T as DTypeIdentifier>::dtype()
    }
}

pub mod topology {
    use super::TopologyT;
    use crate::{
        topology::serial::{SingleTypeEntityTopology, SingleTypeEntityTopologyBorrowed},
        traits::Topology,
    };
    use c_api_tools::concretise_types;

    #[concretise_types(
        gen_type(name = "", replace_with = [""]),
        field(arg = 0, name = "wrap", wrapper = "TopologyT", replace_with = ["SingleTypeEntityTopology", "SingleTypeEntityTopologyBorrowed"]),
    )]
    pub fn topology_sub_entity<T: Topology>(topology: &T, dim: usize, index: usize) -> usize {
        topology.sub_entity(dim, index)
    }

    #[concretise_types(
        gen_type(name = "", replace_with = [""]),
        field(arg = 0, name = "wrap", wrapper = "TopologyT", replace_with = ["SingleTypeEntityTopology", "SingleTypeEntityTopologyBorrowed"]),
    )]
    pub fn topology_sub_entities_size<T: Topology>(topology: &T, dim: usize) -> usize {
        topology.sub_entity_iter(dim).map(|_| 1).sum()
    }

    #[concretise_types(
        gen_type(name = "", replace_with = [""]),
        field(arg = 0, name = "wrap", wrapper = "TopologyT", replace_with = ["SingleTypeEntityTopology", "SingleTypeEntityTopologyBorrowed"]),
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
        field(arg = 0, name = "wrap", wrapper = "TopologyT", replace_with = ["SingleTypeEntityTopology", "SingleTypeEntityTopologyBorrowed"]),
    )]
    pub fn topology_connected_entities_size<T: Topology>(topology: &T, dim: usize) -> usize {
        topology.connected_entity_iter(dim).map(|_| 1).sum()
    }

    #[concretise_types(
        gen_type(name = "", replace_with = [""]),
        field(arg = 0, name = "wrap", wrapper = "TopologyT", replace_with = ["SingleTypeEntityTopology", "SingleTypeEntityTopologyBorrowed"]),
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
        geometry::{SingleElementEntityGeometry, SingleElementEntityGeometryBorrowed},
        traits::{Geometry, Point},
        types::RealScalar,
    };
    use c_api_tools::{concretise_types, DType, DTypeIdentifier};
    use ndelement::ciarlet::CiarletElement;
    use std::ffi::c_void;
    use std::slice::from_raw_parts_mut;

    #[concretise_types(
        gen_type(name = "dtype", replace_with = ["f32", "f64"]),
        field(arg = 0, name = "wrap", wrapper = "GeometryT", replace_with = ["SingleElementEntityGeometry<{{dtype}}, CiarletElement<{{dtype}}>>", "SingleElementEntityGeometryBorrowed<{{dtype}}, CiarletElement<{{dtype}}>>"]),
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
        field(arg = 0, name = "wrap", wrapper = "GeometryT", replace_with = ["SingleElementEntityGeometry<{{dtype}}, CiarletElement<{{dtype}}>>", "SingleElementEntityGeometryBorrowed<{{dtype}}, CiarletElement<{{dtype}}>>"]),
    )]
    pub fn geometry_point_count<G: Geometry>(geometry: &G) -> usize {
        geometry.point_count()
    }

    #[concretise_types(
        gen_type(name = "dtype", replace_with = ["f32", "f64"]),
        field(arg = 0, name = "wrap", wrapper = "GeometryT", replace_with = ["SingleElementEntityGeometry<{{dtype}}, CiarletElement<{{dtype}}>>", "SingleElementEntityGeometryBorrowed<{{dtype}}, CiarletElement<{{dtype}}>>"]),
    )]
    pub fn geometry_degree<G: Geometry>(geometry: &G) -> usize {
        geometry.degree()
    }

    #[concretise_types(
        gen_type(name = "dtype", replace_with = ["f32", "f64"]),
        field(arg = 0, name = "wrap", wrapper = "GeometryT", replace_with = ["SingleElementEntityGeometry<{{dtype}}, CiarletElement<{{dtype}}>>", "SingleElementEntityGeometryBorrowed<{{dtype}}, CiarletElement<{{dtype}}>>"]),
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
        geometry::GeometryMap,
        traits::GeometryMap as GeometryMapTrait,
        types::{Array2D, Array2DBorrowed, RealScalar},
    };
    use c_api_tools::{concretise_types, DType, DTypeIdentifier};
    use std::ffi::c_void;
    use std::slice::from_raw_parts_mut;

    #[concretise_types(
        gen_type(name = "dtype", replace_with = ["f32", "f64"]),
        gen_type(name = "atype", replace_with = ["Array2D<", "Array2DBorrowed<'_, "]),
        field(arg = 0, name = "wrap", wrapper = "GeometryMapT", replace_with = ["GeometryMap<{{dtype}}, {{atype}}{{dtype}}>, {{atype}}usize>>"]),
    )]
    pub fn geometry_map_entity_topology_dimension<GM: GeometryMapTrait>(gmap: &GM) -> usize {
        gmap.entity_topology_dimension()
    }

    #[concretise_types(
        gen_type(name = "dtype", replace_with = ["f32", "f64"]),
        gen_type(name = "atype", replace_with = ["Array2D<", "Array2DBorrowed<'_, "]),
        field(arg = 0, name = "wrap", wrapper = "GeometryMapT", replace_with = ["GeometryMap<{{dtype}}, {{atype}}{{dtype}}>, {{atype}}usize>>"]),
    )]
    pub fn geometry_map_geometry_dimension<GM: GeometryMapTrait>(gmap: &GM) -> usize {
        gmap.geometry_dimension()
    }

    #[concretise_types(
        gen_type(name = "dtype", replace_with = ["f32", "f64"]),
        gen_type(name = "atype", replace_with = ["Array2D<", "Array2DBorrowed<'_, "]),
        field(arg = 0, name = "wrap", wrapper = "GeometryMapT", replace_with = ["GeometryMap<{{dtype}}, {{atype}}{{dtype}}>, {{atype}}usize>>"]),
    )]
    pub fn geometry_map_point_count<GM: GeometryMapTrait>(gmap: &GM) -> usize {
        gmap.point_count()
    }

    #[concretise_types(
        gen_type(name = "dtype", replace_with = ["f32", "f64"]),
        gen_type(name = "atype", replace_with = ["Array2D<", "Array2DBorrowed<'_, "]),
        field(arg = 0, name = "wrap", wrapper = "GeometryMapT", replace_with = ["GeometryMap<{{dtype}}, {{atype}}{{dtype}}>, {{atype}}usize>>"]),
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
        gen_type(name = "atype", replace_with = ["Array2D<", "Array2DBorrowed<'_, "]),
        field(arg = 0, name = "wrap", wrapper = "GeometryMapT", replace_with = ["GeometryMap<{{dtype}}, {{atype}}{{dtype}}>, {{atype}}usize>>"]),
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
        gen_type(name = "atype", replace_with = ["Array2D<", "Array2DBorrowed<'_, "]),
        field(arg = 0, name = "wrap", wrapper = "GeometryMapT", replace_with = ["GeometryMap<{{dtype}}, {{atype}}{{dtype}}>, {{atype}}usize>>"]),
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
        gen_type(name = "atype", replace_with = ["Array2D<", "Array2DBorrowed<'_, "]),
        field(arg = 0, name = "wrap", wrapper = "GeometryMapT", replace_with = ["GeometryMap<{{dtype}}, {{atype}}{{dtype}}>, {{atype}}usize>>"]),
    )]
    pub fn geometry_map_dtype<T: RealScalar + DTypeIdentifier, GM: GeometryMapTrait<T = T>>(
        _gmap: &GM,
    ) -> DType {
        <T as DTypeIdentifier>::dtype()
    }
}
